package org.cchmc.epifluidlab.finaleme.utils;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * PrepareAndSplitInputMatrix
 *
 * Input : CpgMultiMetricsStats.hg19.details.bed.gz (12 columns, header
 * included)
 * Output: out_dir/input_matrix.partXXX.tsv.gz (header included)
 *
 * Assumptions (project-internal, stable):
 * - Input has exactly 12 tab-delimited columns and a header line.
 * - Column order matches CpgMultiMetricsStats output.
 * - methyPrior is either numeric or "NaN" (or empty); NaN/empty is normalized
 * to 0.
 */
public class PrepareAndSplitInputMatrix {

    // 1-based column indices (match CpgMultiMetricsStats output)
    private static final int COL_READNAME_1B = 4; // readName
    private static final int COL_FRAGLEN_1B = 5; // FragLen
    private static final int COL_BASEQ_1B = 9; // baseQ
    private static final int COL_OFFSET_1B = 10; // Offset_frag
    private static final int COL_PRIOR_1B = 12; // methyPrior

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

        PartWriter[] writers = new PartWriter[parts];
        for (int i = 0; i < parts; i++) {
            Path out = outDir.resolve(String.format("input_matrix.part%03d.tsv.gz", i));
            writers[i] = new PartWriter(out);
        }

        // MAX_PER_READ: readName -> emitted lines count
        Map<String, Integer> seen = new HashMap<>(1 << 20);

        long inLines = 0; // data lines read (excluding header)
        long outLines = 0; // data lines written (excluding header)
        long skipped = 0; // filtered or capped

        try (BufferedReader br = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(new FileInputStream(detailsGz.toFile())),
                StandardCharsets.UTF_8))) {

            String header = br.readLine();
            if (header == null)
                throw new IllegalArgumentException("Input is empty: " + detailsGz);

            // propagate header to all parts
            for (PartWriter w : writers)
                w.writeLine(header);

            String line;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '#')
                    continue;
                inLines++;

                // Extract only what we need
                String readName = getColumn(line, COL_READNAME_1B);
                int fragLen = parseInt(getColumn(line, COL_FRAGLEN_1B));
                int baseQ = parseInt(getColumn(line, COL_BASEQ_1B));
                int offset = parseInt(getColumn(line, COL_OFFSET_1B));
                String prior = getColumn(line, COL_PRIOR_1B);

                // Same row-level filters as FinaleMe.processMatrixFile
                if (fragLen <= 30 || fragLen >= 500) {
                    skipped++;
                    continue;
                }
                if (baseQ <= 5) {
                    skipped++;
                    continue;
                }
                if (offset < 0) {
                    skipped++;
                    continue;
                }

                // MAX_PER_READ cap (per fragment/readName)
                int c = seen.getOrDefault(readName, 0) + 1;
                if (c > maxPerRead) {
                    skipped++;
                    continue;
                }
                seen.put(readName, c);

                // Normalize methyPrior: NaN/nan/empty -> 0 (keep everything else as-is)
                if (prior == null || prior.isEmpty() || "NaN".equals(prior) || "nan".equals(prior)) {
                    line = replaceLastColumnWithZero(line);
                }

                int part = (readName.hashCode() & 0x7fffffff) % parts;
                writers[part].writeLine(line);
                outLines++;

                if ((outLines % 1_000_000) == 0) {
                    System.err.println("INFO: written=" + outLines + " read=" + inLines);
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

        System.err.println("INFO: data lines read    = " + inLines);
        System.err.println("INFO: data lines written = " + outLines);
        System.err.println("INFO: skipped lines      = " + skipped);
        System.err.println("INFO: unique readNames   = " + seen.size());
        System.err.println("OK: wrote parts to " + outDir);
    }

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

    /**
     * Extract Nth (1-based) tab-delimited column without splitting whole line.
     */
    private static String getColumn(String line, int col1Based) {
        int target = col1Based - 1;
        int start = 0;
        int col = 0;
        int len = line.length();

        for (int i = 0; i <= len; i++) {
            boolean atEnd = (i == len);
            if (atEnd || line.charAt(i) == '\t') {
                if (col == target)
                    return line.substring(start, i);
                col++;
                start = i + 1;
                if (start > len)
                    return null;
            }
        }
        return null;
    }

    private static int parseInt(String s) {
        // Assumption: internal pipeline produces valid integers here.
        return Integer.parseInt(s);
    }

    /**
     * Replace last column (methyPrior) with "0" (expects 12-column TSV).
     */
    private static String replaceLastColumnWithZero(String line) {
        int lastTab = line.lastIndexOf('\t');
        return line.substring(0, lastTab + 1) + "0";
    }
}
