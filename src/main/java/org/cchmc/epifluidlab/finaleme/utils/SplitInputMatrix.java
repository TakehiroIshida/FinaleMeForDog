package org.cchmc.epifluidlab.finaleme.tools;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class SplitInputMatrix {

    private static final int KEY_COL_1BASED = 4; // readName
    private static final int KEY_COL_0BASED = KEY_COL_1BASED - 1;

    public static void main(String[] args) throws Exception {
        if (args.length < 3) {
            System.err.println(
                    "Usage: java ... SplitInputMatrix <input.tsv.gz> <out_dir> <parts>\n" +
                            "Example: java ... SplitInputMatrix out/input_matrix.tsv.gz out/parts 128");
            System.exit(2);
        }

        Path inputGz = Paths.get(args[0]);
        Path outDir = Paths.get(args[1]);
        int parts = Integer.parseInt(args[2]);

        if (parts <= 1)
            throw new IllegalArgumentException("parts must be >= 2");
        if (!Files.isRegularFile(inputGz))
            throw new FileNotFoundException("Missing input: " + inputGz);

        Files.createDirectories(outDir);

        // 1) open input
        try (BufferedReader br = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(new FileInputStream(inputGz.toFile())),
                StandardCharsets.UTF_8))) {
            String header = br.readLine();
            if (header == null)
                throw new IllegalArgumentException("Input is empty: " + inputGz);

            // 2) prepare temp tsv writers (uncompressed) to avoid too many gzip streams
            BufferedWriter[] writers = new BufferedWriter[parts];
            Path[] tmpTsv = new Path[parts];
            long[] lineCounts = new long[parts]; // data lines only

            for (int i = 0; i < parts; i++) {
                tmpTsv[i] = outDir.resolve(String.format("input_matrix.part%03d.tsv", i));
                writers[i] = Files.newBufferedWriter(tmpTsv[i], StandardCharsets.UTF_8,
                        StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING);
                // write header to each part
                writers[i].write(header);
                writers[i].write('\n');
            }

            // 3) stream lines and dispatch
            String line;
            long totalDataLines = 0;
            while ((line = br.readLine()) != null) {
                if (line.isEmpty())
                    continue;
                if (line.charAt(0) == '#')
                    continue;

                // fast-ish split for key column only
                String readName = getColumn(line, KEY_COL_0BASED);
                if (readName == null || readName.isEmpty()) {
                    // if malformed, skip (or throw)
                    continue;
                }

                int part = (readName.hashCode() & 0x7fffffff) % parts;

                writers[part].write(line);
                writers[part].write('\n');

                lineCounts[part]++;
                totalDataLines++;
                if ((totalDataLines % 1_000_000) == 0) {
                    System.err.println("INFO: dispatched lines=" + totalDataLines);
                }
            }

            // 4) close tsv writers
            for (int i = 0; i < parts; i++) {
                writers[i].flush();
                writers[i].close();
            }

            // 5) gzip each part and remove tmp tsv
            long gzTotalLines = 0;
            for (int i = 0; i < parts; i++) {
                Path gzPath = outDir.resolve(String.format("input_matrix.part%03d.tsv.gz", i));
                gzipFile(tmpTsv[i], gzPath);
                Files.deleteIfExists(tmpTsv[i]);

                System.err.println(String.format("OK: wrote %s (data lines=%d)", gzPath, lineCounts[i]));
                gzTotalLines += lineCounts[i];
            }

            System.err.println("INFO: total data lines=" + totalDataLines);
            System.err.println("INFO: sum(part data lines)=" + gzTotalLines);
            if (totalDataLines != gzTotalLines) {
                throw new IllegalStateException(
                        "Line count mismatch: total=" + totalDataLines + " sumParts=" + gzTotalLines);
            }
        }
    }

    /**
     * Extract Nth (0-based) tab-delimited column without splitting all columns.
     */
    private static String getColumn(String line, int colIndex0) {
        int start = 0;
        int col = 0;
        int len = line.length();

        for (int i = 0; i <= len; i++) {
            boolean atEnd = (i == len);
            if (atEnd || line.charAt(i) == '\t') {
                if (col == colIndex0) {
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

    private static void gzipFile(Path inTsv, Path outGz) throws IOException {
        try (InputStream is = new BufferedInputStream(Files.newInputStream(inTsv));
                OutputStream os = new BufferedOutputStream(new GZIPOutputStream(Files.newOutputStream(outGz,
                        StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)))) {

            byte[] buf = new byte[1024 * 1024];
            int r;
            while ((r = is.read(buf)) >= 0) {
                os.write(buf, 0, r);
            }
        }
    }
}
