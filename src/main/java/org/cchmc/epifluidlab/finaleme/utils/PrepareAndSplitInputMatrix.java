package org.cchmc.epifluidlab.finaleme.utils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class PrepareAndSplitInputMatrix {

    // details / input_matrix の列（1-based）
    private static final int COL_READNAME_1B = 4; // readName
    private static final int COL_FRAGLEN_1B = 5; // FragLen
    private static final int COL_BASEQ_1B = 9; // baseQ
    private static final int COL_OFFSET_1B = 10; // Offset_frag
    private static final int COL_PRIOR_1B = 12; // methyPrior

    // priorの異常WARNを出しすぎない
    private static final long MAX_PRIOR_WARN = 50;

    public static void main(String[] args) throws Exception {
        if (args.length < 4) {
            System.err.println(
                    "Usage:\n" +
                            "  java ... PrepareAndSplitInputMatrix <details.bed.gz> <out_dir> <parts> <max_per_read>\n"
                            +
                            "Example:\n" +
                            "  java ... PrepareAndSplitInputMatrix out/CpgMultiMetricsStats.hg19.details.bed.gz out/parts 128 20");
            System.exit(2);
        }

        Path detailsGz = Paths.get(args[0]);
        Path outDir = Paths.get(args[1]);
        int parts = Integer.parseInt(args[2]);
        int maxPerRead = Integer.parseInt(args[3]);

        if (parts < 2)
            throw new IllegalArgumentException("parts must be >= 2");
        if (maxPerRead < 1)
            throw new IllegalArgumentException("max_per_read must be >= 1");
        if (!Files.isRegularFile(detailsGz))
            throw new FileNotFoundException("Missing input: " + detailsGz);

        Files.createDirectories(outDir);

        // parts個のgzip出力を開く
        PartWriter[] writers = new PartWriter[parts];
        for (int i = 0; i < parts; i++) {
            Path out = outDir.resolve(String.format("input_matrix.part%03d.tsv.gz", i));
            writers[i] = new PartWriter(out);
        }

        // MAX_PER_READ用: readName -> count
        Map<String, Integer> seen = new HashMap<>(1 << 20);

        long totalLines = 0;
        long passedLines = 0;
        long skippedMalformed = 0;
        long skippedFilter = 0;
        long skippedCapped = 0;

        long warnedNonNumericPrior = 0;

        try (BufferedReader br = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(new FileInputStream(detailsGz.toFile())),
                StandardCharsets.UTF_8))) {

            String header = br.readLine();
            if (header == null) {
                throw new IllegalArgumentException("Input is empty: " + detailsGz);
            }

            // ヘッダーを簡易的に確認
            verifyHeaderOrWarn(header);

            // headerを各partに書く
            for (PartWriter w : writers)
                w.writeLine(header);

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty())
                    continue;
                if (line.charAt(0) == '#')
                    continue;

                totalLines++;

                // 12列（タブ11個）チェック
                if (countTabs(line) != 11) {
                    skippedMalformed++;
                    continue;
                }

                // 必要列だけ抜く（readName）
                String readName = getColumn(line, COL_READNAME_1B);
                if (readName == null || readName.isEmpty()) {
                    skippedMalformed++;
                    continue;
                }

                int fragLen = parseIntSafe(getColumn(line, COL_FRAGLEN_1B), Integer.MIN_VALUE);
                int baseQ = parseIntSafe(getColumn(line, COL_BASEQ_1B), Integer.MIN_VALUE);
                int offset = parseIntSafe(getColumn(line, COL_OFFSET_1B), Integer.MIN_VALUE);
                String prior = getColumn(line, COL_PRIOR_1B);

                // FinaleMeと同じフィルタを適用
                if (!(fragLen > 30 && fragLen < 500)) {
                    skippedFilter++;
                    continue;
                }
                if (!(baseQ > 5)) {
                    skippedFilter++;
                    continue;
                }
                if (!(offset >= 0)) {
                    skippedFilter++;
                    continue;
                }

                // MAX_PER_READ
                int c = seen.getOrDefault(readName, 0) + 1;
                if (c > maxPerRead) {
                    skippedCapped++;
                    continue;
                }
                seen.put(readName, c);

                // prior が NaN/empty の場合は 0 に置換
                // 数値でない場合はWARN
                if (prior == null || prior.isEmpty() || "NaN".equals(prior) || "nan".equals(prior)) {
                    line = replaceLastColumnWithZero(line);
                } else if (!looksNumeric(prior)) {
                    // 数値として読めないpriorはWARN
                    if (warnedNonNumericPrior < MAX_PRIOR_WARN) {
                        System.err.println(
                                "WARN: non-numeric methyPrior detected. readName=" + readName + " methyPrior=" + prior);
                    } else if (warnedNonNumericPrior == MAX_PRIOR_WARN) {
                        System.err.println(
                                "WARN: too many non-numeric methyPrior warnings; suppressing further warnings.");
                    }
                    warnedNonNumericPrior++;
                }

                int part = (readName.hashCode() & 0x7fffffff) % parts;
                writers[part].writeLine(line);
                passedLines++;

                if ((passedLines % 1_000_000) == 0) {
                    System.err.println("INFO: wrote lines=" + passedLines + " (seen input=" + totalLines + ")");
                }
            }
        } finally {
            for (PartWriter w : writers) {
                try {
                    w.close();
                } catch (Exception ignored) {
                }
            }
        }

        // サマリ
        System.err.println("INFO: total input lines (non-header) = " + totalLines);
        System.err.println("INFO: passed lines (written)        = " + passedLines);
        System.err.println("INFO: skipped malformed             = " + skippedMalformed);
        System.err.println("INFO: skipped by filters            = " + skippedFilter);
        System.err.println("INFO: skipped by MAX_PER_READ cap   = " + skippedCapped);
        System.err.println("INFO: unique readNames tracked      = " + seen.size());
        System.err.println("INFO: non-numeric methyPrior warns  = " + warnedNonNumericPrior);

        // 出力ファイルの存在とサイズを確認
        for (int i = 0; i < parts; i++) {
            Path out = outDir.resolve(String.format("input_matrix.part%03d.tsv.gz", i));
            if (!Files.isRegularFile(out) || Files.size(out) == 0) {
                throw new IllegalStateException("Missing/empty part: " + out);
            }
        }
        System.err.println("OK: all part files created in " + outDir);
    }

    // ヘッダーの簡易検査用
    private static void verifyHeaderOrWarn(String headerLine) {
        int tabs = countTabs(headerLine);
        if (tabs != 11) {
            System.err.println(
                    "WARN: header does not look like 12 columns (tabCount=" + tabs + "). header=" + headerLine);
            return;
        }

        String[] cols = headerLine.split("\t", -1);
        // 期待列名
        checkHeaderCol(cols, COL_READNAME_1B, "readName");
        checkHeaderCol(cols, COL_FRAGLEN_1B, "FragLen");
        checkHeaderCol(cols, COL_BASEQ_1B, "baseQ");
        checkHeaderCol(cols, COL_OFFSET_1B, "Offset_frag");
        checkHeaderCol(cols, COL_PRIOR_1B, "methyPrior");
    }

    private static void checkHeaderCol(String[] cols, int col1B, String expected) {
        int idx = col1B - 1;
        if (idx < 0 || idx >= cols.length) {
            System.err.println("WARN: header missing column index " + col1B + " (expected=" + expected + ")");
            return;
        }
        String actual = cols[idx];
        if (!expected.equals(actual)) {
            System.err.println(
                    "WARN: header column mismatch at col " + col1B +
                            " expected='" + expected + "' actual='" + actual + "'");
        }
    }

    // gzipファイルに書くBufferedWriterを作るヘルパークラス
    private static final class PartWriter implements Closeable {
        private final BufferedWriter bw;

        PartWriter(Path outGz) throws IOException {
            OutputStream os = new BufferedOutputStream(
                    new GZIPOutputStream(Files.newOutputStream(
                            outGz,
                            StandardOpenOption.CREATE,
                            StandardOpenOption.TRUNCATE_EXISTING)),
                    1 << 20);
            this.bw = new BufferedWriter(new OutputStreamWriter(os, StandardCharsets.UTF_8), 1 << 20);
        }

        void writeLine(String s) throws IOException {
            bw.write(s);
            bw.write('\n');
        }

        @Override
        public void close() throws IOException {
            bw.flush();
            bw.close();
        }
    }

    // 行内のタブ数を確認
    private static int countTabs(String line) {
        int n = 0;
        for (int i = 0, len = line.length(); i < len; i++) {
            if (line.charAt(i) == '\t')
                n++;
        }
        return n;
    }

    // 1-based列を取り出す
    private static String getColumn(String line, int col1Based) {
        int target = col1Based - 1;
        int start = 0;
        int col = 0;
        int len = line.length();

        for (int i = 0; i <= len; i++) {
            boolean atEnd = (i == len);
            if (atEnd || line.charAt(i) == '\t') {
                if (col == target) {
                    return line.substring(start, i);
                }
                col++;
                start = i + 1;
                if (start > len)
                    return null;
            }
        }
        return null;
    }

    private static int parseIntSafe(String s, int fallback) {
        if (s == null)
            return fallback;
        try {
            return Integer.parseInt(s);
        } catch (NumberFormatException e) {
            return fallback;
        }
    }

    // 行の最後の列（methyPrior）を "0" に置換
    private static String replaceLastColumnWithZero(String line) {
        int lastTab = line.lastIndexOf('\t');
        if (lastTab < 0)
            return line;
        return line.substring(0, lastTab + 1) + "0";
    }

    // methyPriorが数値でない場合のWARN抑制用
    private static boolean looksNumeric(String s) {
        if (s == null)
            return false;
        String t = s.trim();
        if (t.isEmpty())
            return false;
        try {
            Double.parseDouble(t);
            return true;
        } catch (NumberFormatException e) {
            return false;
        }
    }
}
