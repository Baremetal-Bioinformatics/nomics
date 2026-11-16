import java.nio.file.*;
import org.broadinstitute.hellbender.utils.io.IOUtils;
/** Exercise GATK 4.6.2.0's archive reader, which stops after one gzip member. */
public final class ExtractArchive {
    public static void main(String[] args) throws Exception {
        Path out = Path.of(args[1]);
        Files.createDirectories(out);
        IOUtils.extractTarGz(Path.of(args[0]), out, true);
        try (var entries = Files.list(out)) {
            long count = entries.filter(Files::isRegularFile).count();
            if (count != 3) throw new AssertionError("expected three members, got " + count);
        }
        System.out.println("PASS GATK archive reader");
    }
}
