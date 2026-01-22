/**
 * CpgMultiMetricsStats.java
 * Feb 27, 2016
 * 5:28:37 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.finaleme.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap; // keyに対応するvalueを格納する辞書型データ構造を生成
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.Logger;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.io.BigWigFileReader;

/**
 *
 */
public class CpgMultiMetricsStats {

	// コマンドラインで使用するオプションの定義（args4jを使用）
	// CpG判定に使う塩基の最小品質スコア
	@Option(name = "-minBaseQ", usage = "minimum base quality score required to check. Default: 5")
	public int minBaseQ = 5;

	// マッピングの最小品質スコア
	@Option(name = "-minMapQ", usage = "minimum mapping quality score required to check. Default: 30")
	public int minMapQ = 30;

	// 最大フラグメント長（異様に長いフラグメントを除外するための閾値）
	@Option(name = "-maxFragLen", usage = "maximum fragment length allowed to check. Default: 500")
	public int maxFragLen = 500;

	// フラグメント端からの最大距離（CpGがフラグメント末端から離れすぎている場合を除去するための閾値）
	@Option(name = "-maxDistToFragEnd", usage = "maximum distant to the end of the fragment allowed to check. in order to be copnsistent with training model. Default: 250")
	public int maxDistToFragEnd = 250;

	// カバレッジ正規化に使用するBAMファイル内の総リード数
	@Option(name = "-totalReadsInBam", usage = "total number of reads used to normalize coverage column. default estimate from bam file by program. Default: -1")
	public long totalReadsInBam = -1;

	// 高カバレッジCpGを除去するための閾値
	@Option(name = "-maxCov", usage = "maximum coverage allowed to check. Default: 250")
	public int maxCov = 250;

	// k-merの最大の長さ
	@Option(name = "-kmerLen", usage = "the K-mer length to check. Default: 0")
	public int kmerLen = 0;

	// 使用するk-merを指定するファイル
	@Option(name = "-kmerString", usage = "the file contain selected K-mer to check. Otherwise, use -kmerLen to automately generate all the k-mer. Default: null")
	public String kmerString = null;

	// CpG周辺の参照配列からk-merを取る範囲
	@Option(name = "-kmerExt", usage = "the +/- region in reference genome to check the k-mer frequency. default is +/- 100bp around CpGs. Default: 100")
	public int kmerExt = 100;

	// 除外する領域を指定するBEDファイル
	@Option(name = "-excludeRegions", usage = "bed files indicated excluded regions. -excludeRegions trackFileName. Default: null")
	public ArrayList<String> excludeRegions = null;

	// 重複領域をチェックするためのBEDファイル
	@Option(name = "-overlapRegions", usage = "bed files to check if regions are overlapped. -overlapRegions trasckName:trackFileName. Default: null")
	public ArrayList<String> overlapRegions = null;

	// CpGから特定の領域までの距離を計算するためのBEDファイル
	@Option(name = "-distantRegions", usage = "bed files to check the distance to these regions. like the distance to TSS. -distantRegions trasckName:trackFileName. Default: null")
	public ArrayList<String> distantRegions = null;

	// bigWigファイルに書かれた数値特徴量をCpG周辺から取得するためのファイル
	@Option(name = "-valueWigs", usage = "bigwig files to check the value in these regions. like the recombination rate of some region.extRegion indicate the +/- regionBp from CpG, -valueWigs trasckName:extRegion:trackFileName. Default: null")
	public ArrayList<String> valueWigs = null;

	// bed.gzに書かれた数値特徴量をCpG周辺から取得するためのファイル
	// trackname: 出力カラム名
	// extRegion: CpGからの拡張領域(bp)
	// trackFileName: tabix index済みのbed.gzファイル
	@Option(name = "-valueBeds", usage = "tabixed bed.gz files to check the value in these regions. like the recombination rate of some region. -valueBeds trasckName:extRegion:trackFileName. Default: null")
	public ArrayList<String> valueBeds = null;

	// WGBSではなくWGSとして扱うモード（Bisulfite判定を無効化）
	@Option(name = "-wgsMode", usage = "used for WGS but not bisulfite space. . Default: false")
	public boolean wgsMode = false;

	// ペアエンドのsecond readを解析から除外
	@Option(name = "-skipSecondEnd", usage = "skip the 2nd end for the statistics. Default: false")
	public boolean skipSecondEnd = false;

	// 隣接CpGまでの距離を特徴量として含めるかどうかのオプション
	@Option(name = "-includeCpgDist", usage = "include the distance between CpGs. Default: false")
	public boolean includeCpgDist = false;

	// 向きが正しい高品質ペアエンドのみを使用するオプション
	@Option(name = "-stringentPaired", usage = "Only use paired end reads that faced to each other. Default: false")
	public boolean stringentPaired = false;

	// 参照配列ではなくフラグメント配列からk-merを計算するオプション
	@Option(name = "-useFragBaseKmer", usage = "add k-mer in fragment. Default: false")
	public boolean useFragBaseKmer = false;

	// フラグメントの向きを考慮してk-merを計算するオプション
	@Option(name = "-useStrandSpecificFragBase", usage = "use k-mer generated by aware the strand of fragment status. Default: false")
	public boolean useStrandSpecificFragBase = false;

	// "chr1"ではなく"1"表記のBAMファイルを使用するオプション
	@Option(name = "-useNoChrPrefixBam", usage = "use bam file with GRch37 instead of hg19 coordinate. Default: false")
	public boolean useNoChrPrefixBam = false;

	@Option(name = "-h", usage = "show option information")
	public boolean help = false;

	// コマンドラインで指定する引数の定義
	@Argument
	private List<String> arguments = new ArrayList<String>();

	// 使用方法の表示
	final private static String USAGE = "CpgMultiMetricsStats [opts] hg19.2bit cpg_list.bed all_cpg.bed wgs.bam cpg_detail.txt.gz";

	private static Logger log = Logger.getLogger(CpgMultiMetricsStats.class);

	private static long startTime = -1;
	private static long points = 0;

	/**
	 * @param args
	 * @throws Exception
	 */
	// Javaプログラムのエントリーポイント
	// CLI引数を受け取り、実処理はdoMainメソッドに委譲
	public static void main(String[] args) throws Exception {
		CpgMultiMetricsStats cmms = new CpgMultiMetricsStats();
		// BasicConfigurator.configure();
		cmms.doMain(args);
	}

	@SuppressWarnings("resource")
	// main()から呼び出され、実際の処理を行うメソッド
	public void doMain(String[] args)
			throws Exception {

		// args4jパーサを生成してコマンドライン引数を解析
		CmdLineParser parser = new CmdLineParser(this);
		// CLI引数のエラー処理
		try {
			if (help || args.length < 5)
				throw new CmdLineException(parser, USAGE, new Throwable());
			// コマンドライン引数を変換・格納
			parser.parseArgument(args);

		} catch (CmdLineException e) {
			System.err.println(e.getMessage());
			parser.printUsage(System.err);
			System.err.println();
			return;
		}

		// 位置引数から入力ファイル・出力ファイルを取得
		String refFile = arguments.get(0); // 参照ゲノム(2bit)
		String cpgListFile = arguments.get(1); // CpGリストファイル(BED)
		String allCpgFile = arguments.get(2); // ゲノム全体のCpGファイル(BED)
		String wgsBamFile = arguments.get(3); // WGS BAMファイル
		String detailFile = arguments.get(4); // HMMに渡すCpG × 特徴量ファイル(gzip)

		// 開始時刻を記録
		initiate();

		// 2bit参照ゲノムファイルの読み込み
		TwoBitParser refParser = new TwoBitParser(new File(refFile));

		// String[] names = p.getSequenceNames();
		// for(int i=0;i<names.length;i++) {
		// p.setCurrentSequence(names[i]);
		// p.printFastaSequence();
		// p.close();
		// }
		// loading exlusion interval file

		// load interval files
		log.info("Processing interval file ... ");
		// BAMファイルの読み込み
		// ValidationStringency.SILENT は BAMの軽微な不整合で止まらないようにする設定
		SamReader wgsReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
				.open(new File(wgsBamFile));
		// SAMSequenceDictionary dictSeq =
		// SAMSequenceDictionaryExtractor.extractDictionary(new File(wgsBamFile));
		// GenomeLocParser glpSeq = new GenomeLocParser(dictSeq);

		// -excludeRegionsで指定されたBEDファイルの読み込み
		// keyに対応する除外区間を格納する辞書型データ構造を生成
		HashMap<String, IntervalTree<Integer>> ignoreLocCollections = null;
		if (excludeRegions != null && !excludeRegions.isEmpty()) {
			log.info("Excluding intervals ... ");
			ignoreLocCollections = new HashMap<String, IntervalTree<Integer>>();

			// 除外領域ファイルを1つずつ処理するループ処理
			for (String excludeRegion : excludeRegions) {
				// BEDファイルを1行ずつ読み込むためのBufferedReaderを生成
				BufferedReader br = new BufferedReader(new FileReader(excludeRegion));
				// 読み込んだ行を格納するための変数
				String line;

				// BEDファイルを1行ずつ読み込むループ処理
				while ((line = br.readLine()) != null) {
					// コメント行はスキップ
					if (line.startsWith("#"))
						continue;
					// タブ区切りで行を分割
					String[] splitLines = line.split("\t");
					// 染色体名・start・endの3つの情報が揃っていない場合はスキップ
					if (splitLines.length < 3) {
						continue;
					}
					// 染色体名を格納
					String chr = splitLines[0];
					// 開始位置を格納
					int start = Integer.parseInt(splitLines[1]);
					// 終了位置を格納
					int end = Integer.parseInt(splitLines[2]);
					IntervalTree<Integer> tree;
					// 既にその染色体名が辞書に存在する場合は対応するIntervalTreeを取得
					if (ignoreLocCollections.containsKey(chr)) {
						tree = ignoreLocCollections.get(chr);
						// 存在しない場合は新しいIntervalTreeを生成
					} else {
						tree = new IntervalTree<Integer>();
					}
					// 除外区間をIntervalTreeに追加
					tree.put(start, end, 1);
					ignoreLocCollections.put(chr, tree);
				}
				br.close();

			}
		}

		// -includeCpgDistオプションがTrueに指定された場合はCpG間の距離を計算（渡されたCpG座標を使用）
		// keyに対応する除外区間を格納する辞書型データ構造を生成（染色体ごとにCpG座標を格納）
		HashMap<String, IntervalTree<String>> allCpgLocCollections = new HashMap<String, IntervalTree<String>>();
		if (includeCpgDist) {
			log.info("Loading all CpG intervals ... ");
			GZIPInputStream gzipInputStream = null;
			BufferedReader br1;
			// .gzファイルの場合はGZIPInputStreamを使用して読み込む
			if (allCpgFile.endsWith(".gz")) {
				gzipInputStream = new GZIPInputStream(new FileInputStream(allCpgFile));
				br1 = new BufferedReader(new InputStreamReader(gzipInputStream));

			} else {
				br1 = new BufferedReader(new FileReader(allCpgFile));
			}

			String line1;

			while ((line1 = br1.readLine()) != null) {
				if (line1.startsWith("#"))
					continue;
				String[] splitLines = line1.split("\t");
				if (splitLines.length < 3) {
					continue;
				}
				String chr = splitLines[0];
				int start = Integer.parseInt(splitLines[1]);
				int end = Integer.parseInt(splitLines[2]);
				IntervalTree<String> tree;

				if (allCpgLocCollections.containsKey(chr)) {
					tree = allCpgLocCollections.get(chr);
				} else {
					tree = new IntervalTree<String>();
				}
				// strandはゲノム配列に対する向きを示す
				String strand = ".";
				// strandがある場合は取得
				if (splitLines.length >= 6) {
					if (splitLines[5].equalsIgnoreCase("-")) {
						strand = "-";
					} else if (splitLines[5].equalsIgnoreCase("+")) {
						strand = "+";
					}
					// strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
				}
				tree.put(start, end, strand);
				allCpgLocCollections.put(chr, tree);
			}
			if (allCpgFile.endsWith(".gz")) {
				gzipInputStream.close();
			}
			br1.close();
		}

		// -overlapRegionsで指定されたBEDファイルの読み込み
		// CpGがBED区間と重複するかどうかを追加カラムとして出力するために使用(0 or 1)
		HashMap<String, HashMap<String, IntervalTree<Integer>>> overlapLocStringCollections = null;
		LinkedHashSet<String> overlapLocString = new LinkedHashSet<String>();
		if (overlapRegions != null && !overlapRegions.isEmpty()) {
			log.info("Overalpped intervals ... ");
			overlapLocStringCollections = new HashMap<String, HashMap<String, IntervalTree<Integer>>>();
			for (String overlapRegionString : overlapRegions) {

				String[] splitStrings = overlapRegionString.split(":");
				if (splitStrings.length < 2) {
					throw new IllegalArgumentException("need to provide trackname:trackFile for overlapRegions");
				}
				HashMap<String, IntervalTree<Integer>> overlapLocCollections = new HashMap<String, IntervalTree<Integer>>();
				String overlapRegionName = splitStrings[0];
				String overlapRegion = splitStrings[1];
				overlapLocString.add(overlapRegionName);
				BufferedReader br = new BufferedReader(new FileReader(overlapRegion));
				String line;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("#"))
						continue;
					String[] splitLines = line.split("\t");
					if (splitLines.length < 3) {
						continue;
					}
					String chr = splitLines[0];
					int start = Integer.parseInt(splitLines[1]);
					int end = Integer.parseInt(splitLines[2]);
					IntervalTree<Integer> tree;
					if (overlapLocCollections.containsKey(chr)) {
						tree = overlapLocCollections.get(chr);
					} else {
						tree = new IntervalTree<Integer>();
					}
					tree.put(start, end, 1);
					overlapLocCollections.put(chr, tree);

				}
				br.close();
				overlapLocStringCollections.put(overlapRegionName, overlapLocCollections);

			}
		}

		// -distantRegionsで指定されたBEDファイルの読み込み
		// CpGからBED区間までの距離を追加カラムとして出力するために使用
		HashMap<String, HashMap<String, IntervalTree<String>>> distantLocStringCollections = null;
		LinkedHashSet<String> distantLocString = new LinkedHashSet<String>();
		if (distantRegions != null && !distantRegions.isEmpty()) {
			log.info(" Intervals used to calculate distances... ");
			distantLocStringCollections = new HashMap<String, HashMap<String, IntervalTree<String>>>();
			for (String distantRegionString : distantRegions) {

				String[] splitStrings = distantRegionString.split(":");
				if (splitStrings.length < 2) {
					throw new IllegalArgumentException("need to provide trackname:trackFile for distantRegions");
				}
				HashMap<String, IntervalTree<String>> distantLocCollections = new HashMap<String, IntervalTree<String>>();
				String distantRegionName = splitStrings[0];
				String distantRegion = splitStrings[1];
				distantLocString.add(distantRegionName);

				BufferedReader br = new BufferedReader(new FileReader(distantRegion));
				String line;
				while ((line = br.readLine()) != null) {
					if (line.startsWith("#"))
						continue;
					String[] splitLines = line.split("\t");
					if (splitLines.length < 3) {
						continue;
					} else {
						String chr = splitLines[0];
						int start = Integer.parseInt(splitLines[1]);
						int end = Integer.parseInt(splitLines[2]);
						IntervalTree<String> tree;
						if (distantLocCollections.containsKey(chr)) {
							tree = distantLocCollections.get(chr);
						} else {
							tree = new IntervalTree<String>();
						}
						String strand = ".";
						if (splitLines.length >= 6) {
							if (splitLines[5].equalsIgnoreCase("-")) {
								strand = "-";
							} else if (splitLines[5].equalsIgnoreCase("+")) {
								strand = "+";
							}
							// strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
						}
						tree.put(start, end, strand);
						distantLocCollections.put(chr, tree);
					}

				}
				br.close();
				distantLocStringCollections.put(distantRegionName, distantLocCollections);

			}
		}

		HashMap<String, Pair<Integer, TabixFeatureReader<BEDFeature, ?>>> valueBedReaders = null;
		LinkedHashSet<String> valueBedLocString = new LinkedHashSet<String>();
		if (valueBeds != null) {
			log.info("Loading value interval bed file ... ");
			valueBedReaders = new HashMap<String, Pair<Integer, TabixFeatureReader<BEDFeature, ?>>>();
			for (String valueBedString : valueBeds) {
				String[] splitStrings = valueBedString.split(":");
				if (splitStrings.length < 3) {
					throw new IllegalArgumentException("need to provide trackname:extRegion:trackFile for valueBeds");
				}
				String valueBedName = splitStrings[0];
				int valueBedExt = Integer.parseInt(splitStrings[1]);
				String valueRegion = splitStrings[2];
				valueBedLocString.add(valueBedName);
				valueBedReaders.put(valueBedName, new Pair<Integer, TabixFeatureReader<BEDFeature, ?>>(valueBedExt,
						new TabixFeatureReader(valueRegion, new BEDCodec())));
			}

		}

		HashMap<String, Pair<Integer, BigWigFileReader>> valueWigReaders = null;
		LinkedHashSet<String> valueWigLocString = new LinkedHashSet<String>();
		if (valueWigs != null) {
			log.info("Loading value interval big wig file ... ");
			valueWigReaders = new HashMap<String, Pair<Integer, BigWigFileReader>>();
			for (String valueWigString : valueWigs) {
				String[] splitStrings = valueWigString.split(":");
				if (splitStrings.length < 3) {
					throw new IllegalArgumentException("need to provide trackname:extRegion:trackFile for valueWigs");
				}
				String valueWigName = splitStrings[0];
				int valueWigExt = Integer.parseInt(splitStrings[1]);
				String valueRegion = splitStrings[2];
				valueWigLocString.add(valueWigName);
				valueWigReaders.put(valueWigName, new Pair<Integer, BigWigFileReader>(valueWigExt,
						new BigWigFileReader(new File(valueRegion).toPath())));
			}

		}

		// k-merの準備
		// -kmerStringオプションで指定されたファイルからk-merを読み込むか、kmerLenオプションで指定された長さまでの全k-merを自動生成
		// 最終的に列として出力するk-merの集合
		LinkedHashSet<String> kmerCollections = new LinkedHashSet<String>(); // 重複なしで順序を保持する集合
		// k-merファイルが指定されていればそれを読み込む
		if (kmerString != null) {
			log.info("Loading selected K-mer file ... ");
			BufferedReader br = new BufferedReader(new FileReader(kmerString));
			String line;
			while ((line = br.readLine()) != null) {
				if (line.startsWith("#"))
					continue;
				kmerCollections.add(line);

			}
			br.close();
			// k-merファイルが指定されていなければkmerLenまでの全k-merを自動生成
		} else {
			log.info("Automate generate all k-mer until length " + kmerLen);
			for (int i = 2; i <= kmerLen; i++) {
				for (byte[] kmer : SequenceUtil.generateAllKmers(i)) {
					kmerCollections.add(new String(kmer));

				}
			}

		}

		// 出力ファイルの追加カラム名を準備
		String header = "";
		if (overlapLocString.size() > 0) {
			for (String key : overlapLocString) {
				header = header + "\t" + key;
			}
		}
		if (distantLocString.size() > 0) {
			for (String key : distantLocString) {
				header = header + "\t" + key;
			}
		}
		if (valueBedLocString.size() > 0) {
			for (String key : valueBedLocString) {
				header = header + "\t" + key;
			}
		}
		if (valueWigLocString.size() > 0) {
			for (String key : valueWigLocString) {
				header = header + "\t" + key;
			}
		}
		if (kmerCollections.size() > 0) {
			if (!useFragBaseKmer) {
				for (String key : kmerCollections) {
					header = header + "\t" + key;
				}
			} else {
				for (String key : kmerCollections) {
					header = header + "\t" + key + "_frag";
				}
			}

		}

		// 解析対象CpGリストの読み込み
		log.info("Loading CpG interval file ... ");
		HashMap<String, IntervalTree<String>> cpgCollections = new HashMap<String, IntervalTree<String>>();
		GZIPInputStream gzipInputStream1 = null;
		BufferedReader br;
		if (allCpgFile.endsWith(".gz")) {
			gzipInputStream1 = new GZIPInputStream(new FileInputStream(cpgListFile));
			br = new BufferedReader(new InputStreamReader(gzipInputStream1));

		} else {
			br = new BufferedReader(new FileReader(cpgListFile));
		}

		String line;

		while ((line = br.readLine()) != null) {
			if (line.startsWith("#"))
				continue;
			String[] splitLines = line.split("\t");
			if (splitLines.length < 3) {
				continue;
			}
			String chr = splitLines[0];
			int start = Integer.parseInt(splitLines[1]);
			int end = Integer.parseInt(splitLines[2]);
			IntervalTree<String> tree;

			if (cpgCollections.containsKey(chr)) {
				tree = cpgCollections.get(chr);
			} else {
				tree = new IntervalTree<String>();
			}
			String strand = ".";
			if (splitLines.length >= 6) {
				if (splitLines[5].equalsIgnoreCase("-")) {
					strand = "-";
				} else if (splitLines[5].equalsIgnoreCase("+")) {
					strand = "+";
				}
				// strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
			}
			tree.put(start, end, strand);
			cpgCollections.put(chr, tree);
		}
		if (cpgListFile.endsWith(".gz")) {
			gzipInputStream1.close();
		}
		br.close();

		// 正規化用の総リード数を計算
		// -totalReadsInBamオプションが指定されていればそれを使用、そうでなければBAMファイルから計算
		double readsNumTotal = 0;
		if (totalReadsInBam > 0) {
			log.info("Get total reads number used for scaling from input option -totalReadsInBam ... ");
			readsNumTotal = totalReadsInBam;
		} else {
			log.info("Get total reads number used for scaling from bam file... ");
			SAMRecordIterator wgsIt = wgsReader.iterator();

			while (wgsIt.hasNext()) {
				SAMRecord r = wgsIt.next();
				if (failFlagFilter(r)) {
					continue;
				} else {
					if (stringentPaired && !CcInferenceUtils.passReadPairOrientation(r)) {
						continue;
					}

				}
				readsNumTotal++;
			}
			wgsIt.close();
		}

		log.info((long) readsNumTotal + " reads in total ...");
		readsNumTotal = readsNumTotal / 1000000;
		log.info("Output value for each CpG in each DNA fragment ... ");
		FileOutputStream output = new FileOutputStream(detailFile);
		// 出力ファイルを開いてヘッダーを書き込む
		OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");

		// basic
		writer.write(
				"chr\tstart\tend\treadName\tFragLen\tFrag_strand\tmethy_stat\tNorm_Frag_cov\tbaseQ\tOffset_frag\tDist_frag_end");
		if (includeCpgDist) {
			writer.write("\tdist_nearest_CpG");
		}
		writer.write(header + "\n");
		// fragment k-mer

		//

		long i = 0;
		String prevChr = "";
		// cpgCollectionsは解析対象CpGのIntervalTreeをchrごとに持っている（辞書型データ構造）
		for (String chr : cpgCollections.keySet()) {
			if (chr.equalsIgnoreCase("chrM")) {
				continue;

			}
			if (prevChr.equalsIgnoreCase("")) {
				prevChr = chr;
				refParser.setCurrentSequence(chr);
			}
			String bamChr = chr;
			if (useNoChrPrefixBam) {

				Pattern replace = Pattern.compile("^chr");
				Matcher matcher1 = replace.matcher(bamChr);
				bamChr = matcher1.replaceAll("");
			}
			IntervalTree<String> cpgChrCollections = cpgCollections.get(chr);
			Iterator<Node<String>> cpgIterator = cpgChrCollections.iterator();
			while (cpgIterator.hasNext()) {
				Node<String> cpg = cpgIterator.next();
				int start = cpg.getStart();
				int end = cpg.getEnd();
				int fragMostLeft = start + 1;
				int fragMostRight = end;

				SAMRecordIterator wgsIt = wgsReader.queryOverlapping(bamChr, start + 1, end);// start 1-based, inclusive
																								// start of interval of
																								// interest. Zero
																								// implies start of the
																								// reference sequence.
				HashMap<String, SAMRecord> countedReads = new HashMap<String, SAMRecord>();
				// log.info("testincpg" + i + "\t" + cpgCollections.size() + "\t" +
				// wgsIt.hasNext() + "\t" + chr + "\t" + start + "\t" + end);
				int readNumber = 0;
				// log.info(r.getReadName() + "\t" + failFlagFilter(r));
				while (wgsIt.hasNext()) {
					SAMRecord r = wgsIt.next();
					// log.info(r.getReadName() + "\t" + failFlagFilter(r));
					if (failFlagFilter(r)) {
						continue;
					} else {
						if (stringentPaired && !CcInferenceUtils.passReadPairOrientation(r)) {
							continue;
						}
						readNumber++;
						boolean negStrand = r.getReadNegativeStrandFlag();
						boolean secondEnd = r.getReadPairedFlag() && r.getSecondOfPairFlag();
						if (secondEnd) {
							negStrand = !negStrand;
						}

						int bisulfitePos = 0;
						if (!wgsMode) {
							if (r.getTransientAttribute("BS") != null) { // if the reads had been processed before
								bisulfitePos = Integer.parseInt((String) r.getTransientAttribute("BS"));
							} else {
								bisulfitePos = CcInferenceUtils.bisulfiteIncompleteReads(r);
								r.setTransientAttribute("BS", bisulfitePos);

							}
						}

						// log.info("testin" + readsNumTotal + "\t" + bisulfitePos);
						int offSet = r.getReadPositionAtReferencePosition(end) - 1;
						if (bisulfitePos < 0) {
							continue;
						} else if (bisulfitePos > 0) {
							if ((!negStrand && offSet < bisulfitePos) || (negStrand && offSet >= bisulfitePos)) {
								continue;
							}
						}
						if (offSet < 0) {
							// log.info(r.getReadName() + "\t" + offSet + "\t" + end + "\t" +
							// r.getAlignmentStart() + "\t" + r.getAlignmentEnd());
							continue;
						}

						if (r.getAlignmentStart() < fragMostLeft) {
							fragMostLeft = r.getAlignmentStart();
						}
						if (r.getMateAlignmentStart() < fragMostLeft) {
							fragMostLeft = r.getMateAlignmentStart();
						}

						if (r.getAlignmentEnd() > fragMostRight) {
							fragMostRight = r.getAlignmentEnd();
						}
						int mateEnd = CcInferenceUtils.getMateAlignmentEndByMateCigar(r);
						if (mateEnd > fragMostRight) {
							fragMostRight = mateEnd;
						}

						String readName = r.getReadName();
						if (countedReads.containsKey(readName)) {// to filter overlapped fragments, which affect a lot
																	// in cfDNA
							SAMRecord prev = countedReads.get(readName);
							int offSetPrev = prev.getReadPositionAtReferencePosition(end) - 1;
							if (offSet < r.getBaseQualities().length && offSetPrev < prev.getBaseQualities().length) {
								byte baseQ = r.getBaseQualities()[offSet];
								byte base = CcInferenceUtils.toUpperCase(r.getReadBases()[offSet]);

								byte baseQPrev = prev.getBaseQualities()[offSetPrev];
								byte basePrev = CcInferenceUtils.toUpperCase(prev.getReadBases()[offSetPrev]);

								if (!BaseUtils.basesAreEqual(base, basePrev)) {
									if (baseQ > baseQPrev) {
										countedReads.put(readName, r);
									} else if (baseQ < baseQPrev) {

									} else {
										if (!secondEnd) {
											countedReads.put(readName, r);
										}
									}
								}
							}

						} else {
							countedReads.put(readName, r);
						}
					}

				}
				wgsIt.close();
				// if(end==152011200){
				// log.info(countedReads.size() + " reads finally used ...");
				// }
				// log.info(countedReads.size() + " reads finally used ...");
				if (readNumber >= maxCov || countedReads.size() == 0) {
					continue;
				}
				// log.info("test" + readsNumTotal + "\t" + countedReads.size());
				double normalizedFragCov = (double) readNumber / readsNumTotal;
				// System.err.println(normalizedFragCov + "\t" + countedReads.size() + "\t" +
				// readsNumTotal);

				// System.err.println(normalizedFragCov + "\t" + countedReads.size() + "\t" +
				// readsNumTotal);
				if (!chr.equalsIgnoreCase(prevChr)) {
					refParser.close();
					refParser.setCurrentSequence(chr);
					prevChr = chr;
				}
				byte[] refBasesExt = CcInferenceUtils
						.toUpperCase(refParser.loadFragment(end - 1 - kmerExt, kmerExt * 2 + 1).getBytes());
				byte refBase = refBasesExt[kmerExt];

				// if(end==152011200){
				// log.info(new String(refBasesExt) + "\t" + (end-1-kmerExt) + "\t" +
				// (kmerExt*2+1) + "\t" + kmerExt + "\t" + (char)refBasesExt[kmerExt]);

				// }

				HashMap<String, Double> kmerMapsRef = new HashMap<String, Double>();
				if (!useFragBaseKmer) {
					for (int j = 2; j <= kmerLen; j++) {
						kmerMapsRef.putAll(CcInferenceUtils.kmerFreqSearch(refBasesExt, j));

					}
				}

				// nearest cpg's distance in reference genome

				IntervalTree<String> cpgLocCollections = null;
				double nearestCpg = Double.NaN;
				if (includeCpgDist) {
					cpgLocCollections = allCpgLocCollections.get(chr);
					Iterator<Node<String>> upstreamCpgIt = null;
					Iterator<Node<String>> downstreamCpgIt = null;
					if (BaseUtils.basesAreEqual(refBase, BaseUtilsMore.C)) {
						upstreamCpgIt = cpgLocCollections.reverseIterator(start - 1, end - 1);
						downstreamCpgIt = cpgLocCollections.iterator(start + 2, end + 2);

					} else if (BaseUtils.basesAreEqual(refBase, BaseUtilsMore.G)) {
						upstreamCpgIt = cpgLocCollections.reverseIterator(start - 2, end - 2);
						downstreamCpgIt = cpgLocCollections.iterator(start + 1, end + 1);
					} else {
						continue;
					}

					if (upstreamCpgIt == null || !upstreamCpgIt.hasNext()) {
						IntervalTree.Node<String> downstream = downstreamCpgIt.next();
						nearestCpg = CcInferenceUtils.intervalDistance(downstream, cpg);
						// System.err.println(downstream.toString());
					} else if (downstreamCpgIt == null || !downstreamCpgIt.hasNext()) {
						IntervalTree.Node<String> upstream = upstreamCpgIt.next();
						nearestCpg = CcInferenceUtils.intervalDistance(upstream, cpg);
						// System.err.println(upstream.toString());
					} else {
						IntervalTree.Node<String> upstream = upstreamCpgIt.next();
						IntervalTree.Node<String> downstream = downstreamCpgIt.next();
						// System.err.println(upstream.toString());
						// System.err.println(downstream.toString());

						int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
						int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
						if (Math.abs(dist1) < Math.abs(dist2)) {
							nearestCpg = dist1;
						} else {
							nearestCpg = dist2;
						}
					}
				}

				// overlap with feature in reference genome

				HashMap<String, Integer> overlapStatCollections = new HashMap<String, Integer>();
				if (overlapRegions != null && !overlapRegions.isEmpty()) {
					for (String key : overlapLocStringCollections.keySet()) {
						HashMap<String, IntervalTree<Integer>> tmp = overlapLocStringCollections.get(key);
						if (tmp.containsKey(chr)) {
							if (tmp.get(chr).minOverlapper(start, end) == null) {
								overlapStatCollections.put(key, 0);
							} else {
								overlapStatCollections.put(key, 1);
							}

						} else {
							overlapStatCollections.put(key, 0);
						}

					}
				}

				// distance with feature in reference genome
				HashMap<String, Integer> distStatCollections = new HashMap<String, Integer>();
				if (distantRegions != null && !distantRegions.isEmpty()) {
					for (String key : distantLocStringCollections.keySet()) {
						IntervalTree<String> locCollections = distantLocStringCollections.get(key).get(chr);

						// IntervalTree.Node<String> upstream = locCollections.max(start, end);
						// IntervalTree.Node<String> downstream = locCollections.min(start, end);
						int distanceNearest = Integer.MAX_VALUE;
						if (locCollections != null && locCollections.size() > 0) {
							Iterator<Node<String>> upstreamIt = locCollections.reverseIterator(start, end);
							Iterator<Node<String>> downstreamIt = locCollections.iterator(start, end);
							if (!upstreamIt.hasNext()) {
								IntervalTree.Node<String> downstream = locCollections.min(start, end);
								distanceNearest = CcInferenceUtils.intervalDistance(downstream, cpg);
								// System.err.println(downstream.toString());
							} else if (!downstreamIt.hasNext()) {
								IntervalTree.Node<String> upstream = locCollections.max(start, end);
								distanceNearest = CcInferenceUtils.intervalDistance(upstream, cpg);
								// System.err.println(upstream.toString());
							} else {
								IntervalTree.Node<String> upstream = locCollections.max(start, end);
								IntervalTree.Node<String> downstream = locCollections.min(start, end);
								// System.err.println(upstream.toString());
								// System.err.println(downstream.toString());

								int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
								int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
								if (Math.abs(dist1) < Math.abs(dist2)) {
									distanceNearest = dist1;
								} else {
									distanceNearest = dist2;
								}
							}
						}

						distStatCollections.put(key, distanceNearest);
					}
				}

				// value in bed file
				// key: trackname
				// value: Integer = 拡張領域(bp), TabixFeatureReader<BEDFeature, ?> = 高速区間検索用リーダー
				HashMap<String, Double> valBedStatCollections = new HashMap<String, Double>();
				if (valueBedReaders != null) {
					for (String key : valueBedReaders.keySet()) {
						int range = valueBedReaders.get(key).getFirst();
						boolean mean0 = range < 0 ? true : false;
						if (mean0) {
							range = 0 - range;
						}
						TabixFeatureReader<BEDFeature, ?> bedReader = valueBedReaders.get(key).getSecond();
						// CpGの位置を中心にrange bp拡張した区間のBED特徴量を取得
						CloseableTribbleIterator<BEDFeature> featureIt = bedReader.query(chr,
								(start - range < 0 ? 1 : start - range + 1), end + range);
						DescriptiveStatistics statFeature = new DescriptiveStatistics();
						while (featureIt.hasNext()) {
							BEDFeature term = featureIt.next();
							// BEDのscore値を統計量に追加
							if (!Double.isNaN(term.getScore())) {
								statFeature.addValue(term.getScore());
							}

						}
						featureIt.close();
						if (statFeature.getN() > 0) {
							if (mean0) {
								valBedStatCollections.put(key,
										(double) statFeature.getSum() / (double) (range * 2 + 1));
							} else {
								valBedStatCollections.put(key, statFeature.getMean());
							}

						} else {
							valBedStatCollections.put(key, Double.NaN);
						}

					}
				}

				// value in wig file
				// key: trackname
				// value: Integer = 拡張領域(bp), BigWigFileReader = 高速区間検索用リーダー
				HashMap<String, Double> valWigStatCollections = new HashMap<String, Double>();
				if (valueWigReaders != null) {
					for (String key : valueWigReaders.keySet()) {
						int range = valueWigReaders.get(key).getFirst();
						if (range < 0) {
							BigWigFileReader wigReader = valueWigReaders.get(key).getSecond();
							range = 0 - range;
							// queryStatsでmean/sum/Nを取得
							SummaryStatistics statFeature = wigReader.queryStats(chr,
									(start - range < 0 ? 1 : start - range), end + range);

							if (statFeature.getN() > 0) {
								valWigStatCollections.put(key,
										(double) statFeature.getSum() / (double) (range * 2 + 1));
							} else {
								valWigStatCollections.put(key, Double.NaN);
							}

						} else {
							BigWigFileReader wigReader = valueWigReaders.get(key).getSecond();
							SummaryStatistics statFeature = wigReader.queryStats(chr, start - range, end + range);
							if (statFeature.getN() > 0) {
								valWigStatCollections.put(key, statFeature.getMean());
							} else {
								valWigStatCollections.put(key, Double.NaN);
							}
						}

					}
				}

				for (String readName : countedReads.keySet()) {
					SAMRecord r = countedReads.get(readName);
					boolean negStrand = r.getReadNegativeStrandFlag();
					boolean secondEnd = r.getReadPairedFlag() && r.getSecondOfPairFlag();
					if (secondEnd) {
						negStrand = !negStrand;
					}
					int offSet = r.getReadPositionAtReferencePosition(end) - 1;
					if (offSet < 0) { // even it is within the reference interval, but it might be Deletion in the
										// reads
						continue;
					}
					char methyStat = '.';
					byte[] bases = r.getReadBases();

					// System.err.println(readName + "\t" + offSet + "\t" + new String(bases) + "\t"
					// + negStrand + "\t" + secondEnd + "\t" + end + "\t" + r.getAlignmentStart() +
					// "\t" + r.getAlignmentEnd());
					byte base = bases[offSet];
					byte[] baseQs = r.getBaseQualities();
					byte baseQ = baseQs[offSet];
					if (baseQ <= minBaseQ) {
						continue;
					}
					// if(end==152011200){
					// log.info(readName + "\t" + refBase);
					// }
					if (negStrand) {
						// if(end==152011200){
						// log.info(refBase + "\t" + base + "\t" + BaseUtils.basesAreEqual(refBase,
						// BaseUtilsMore.G) + "\t" + BaseUtils.basesAreEqual(base, BaseUtilsMore.G) +
						// "\t" + BaseUtils.basesAreEqual(base, BaseUtilsMore.A));
						// }
						if (BaseUtils.basesAreEqual(refBase, BaseUtilsMore.G)) {
							if (BaseUtils.basesAreEqual(base, BaseUtilsMore.G)) {
								methyStat = 'm';
							} else if (BaseUtils.basesAreEqual(base, BaseUtilsMore.A)) {
								methyStat = 'u';
							} else {
								continue;
							}
						} else {
							continue;
						}
					} else {
						if (end == 152011200) {
							// log.info(refBase + "\t" + base + "\t" + BaseUtils.basesAreEqual(refBase,
							// BaseUtilsMore.C) + "\t" + BaseUtils.basesAreEqual(base, BaseUtilsMore.C) +
							// "\t" + BaseUtils.basesAreEqual(base, BaseUtilsMore.T));
						}
						if (BaseUtils.basesAreEqual(refBase, BaseUtilsMore.C)) {
							if (BaseUtils.basesAreEqual(base, BaseUtilsMore.C)) {
								methyStat = 'm';
							} else if (BaseUtils.basesAreEqual(base, BaseUtilsMore.T)) {
								methyStat = 'u';
							} else {
								continue;
							}
						} else {
							continue;
						}
					}

					// if(end==152011200){
					// log.info("after " + readName);
					// }

					int fragLen = Math.abs(r.getInferredInsertSize());
					if (fragLen > maxFragLen) {
						continue;
					}
					int cpgOffset = CcInferenceUtils.getFragOffsetFromReadsOffset(r, offSet);
					// int distToFragEnd = Math.min((fragLen-cpgOffset), cpgOffset);
					int distToFragEnd = CcInferenceUtils.getDistFragEndFromReadsOffset(r, offSet);
					if (distToFragEnd > maxDistToFragEnd) {
						continue;
					}
					// if(negStrand){
					// offSet = fragLen-offSet;
					// }

					char fragStrand = negStrand ? '-' : '+';

					// get k-mer for the fragment
					int fragStart = Math.min(r.getAlignmentStart(), r.getAlignmentStart() + r.getInferredInsertSize());
					int fragEnd = Math.max(r.getAlignmentStart(), r.getAlignmentStart() + r.getInferredInsertSize());
					if (r.getInferredInsertSize() == 0) {
						continue;
					}

					// byte[] refBasesFragAll = refParser.loadFragment(fragMostLeft,
					// fragMostRight-fragMostLeft+1).getBytes();
					// byte[] refBasesFrag = new byte[fragLen]; //TODO check if it is correct
					// here!!!!
					// System.err.println(r.getAlignmentStart() + "\t" + r.getMateAlignmentStart() +
					// "\t" + fragLen + "\t" + refBasesFragAll.length + "\t" + fragStart + "\t" +
					// fragMostLeft + "\t" + r.getReadName());
					// for(int j = fragStart-fragMostLeft, index = 0; index < fragLen; j++,
					// index++){
					// refBasesFrag[index] = refBasesFragAll[j];
					// }
					byte[] refBasesFrag = refParser.loadFragment(fragMostLeft, fragMostRight - fragMostLeft + 1)
							.getBytes();
					if (negStrand && useStrandSpecificFragBase) {
						SequenceUtil.reverseComplement(refBasesFrag);
					}
					HashMap<String, Double> kmerMapsFrag = new HashMap<String, Double>();
					if (useFragBaseKmer) {
						for (int j = 2; j <= kmerLen; j++) {
							kmerMapsFrag.putAll(CcInferenceUtils.kmerFreqSearch(refBasesFrag, j));
						}
					}

					// System.err.println(CcInferenceUtils.getFragOffsetFromReadsOffset(r, offSet));
					writer.write(chr + "\t" + start + "\t" + end + "\t" + readName + "\t" + fragLen + "\t" + fragStrand
							+ "\t" + methyStat + "\t" + String.format("%.6f", normalizedFragCov)
							+ "\t" + (int) baseQ + "\t" + cpgOffset + "\t" + distToFragEnd);
					if (includeCpgDist) {
						writer.write("\t" + nearestCpg);
					}
					// overlap regions
					if (overlapStatCollections.size() > 0) {
						for (String key : overlapLocString) {
							writer.write("\t" + overlapStatCollections.get(key));
						}
					}

					// distant regions
					if (distStatCollections.size() > 0) {
						for (String key : distantLocString) {
							writer.write("\t" + distStatCollections.get(key));
						}
					}

					// valBed regions
					if (valBedStatCollections.size() > 0) {
						for (String key : valueBedLocString) {
							writer.write("\t" + String.format("%.3f", valBedStatCollections.get(key)));
						}
					}

					// valWig regions
					if (valWigStatCollections.size() > 0) {
						for (String key : valueWigLocString) {
							writer.write("\t" + String.format("%.3f", valWigStatCollections.get(key)));
						}
					}
					// k-mer in reference genome
					// if(kmerMapsRef.size()>0 && useRefSeqBaseKmer){
					if (kmerMapsRef.size() > 0) {
						for (String key : kmerCollections) {
							writer.write("\t" + String.format("%.3f", kmerMapsRef.get(key)));
						}
					}

					// k-mer in fragment
					// if(kmerMapsFrag.size()>0 && !useRefSeqBaseKmer){
					if (kmerMapsFrag.size() > 0) {
						for (String key : kmerCollections) {
							writer.write("\t" + String.format("%.3f", kmerMapsFrag.get(key)));
						}
					}

					writer.write("\n");
					points++;

				}
				i++;
				if (i % 1000 == 0) {
					log.info("Processing Cpg " + i + " ...");
					writer.flush();
				}
			}
		}
		writer.close();
		output.close();

		wgsReader.close();
		refParser.closeParser();
		;

		if (valueBedReaders != null) {
			for (String key : valueBedReaders.keySet()) {
				valueBedReaders.get(key).getSecond().close();
			}
		}
		if (valueWigReaders != null) {
			for (String key : valueWigReaders.keySet()) {
				valueWigReaders.get(key).getSecond().close();
			}
		}

		// 終了時刻を記録
		finish();

	}

	private boolean failFlagFilter(SAMRecord r) {
		return r.getReadUnmappedFlag() || r.getNotPrimaryAlignmentFlag() || r.getMappingQuality() < minMapQ
				|| r.getReadFailsVendorQualityCheckFlag() || r.getDuplicateReadFlag() || !r.getReadPairedFlag()
				|| !r.getProperPairFlag()
				|| (skipSecondEnd && r.getReadPairedFlag() && r.getSecondOfPairFlag());
	}

	// 処理開始時刻を記録
	private void initiate() {
		startTime = System.currentTimeMillis();
	}

	// 処理終了時刻を記録して処理時間をログに出力
	private void finish() {
		long endTime = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime / 60;
		double totalTimeHours = totalTime / 3600;

		log.info("Counted " + points + " data points in total");
		log.info("CpgMultiMetricsStats's running time is: " + String.format("%.2f", totalTime) + " secs, "
				+ String.format("%.2f", totalTimeMins) + " mins, " + String.format("%.2f", totalTimeHours) + " hours");
	}

}
