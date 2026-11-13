import java.io.*;
import java.nio.file.*;
import org.apache.commons.compress.archivers.tar.*;
import org.apache.commons.compress.compressors.gzip.*;
public class ArchiveOracle {
 public static void main(String[] args)throws Exception {
  // Same stream classes and constructor chain as GATK IOUtils.extractTarGz.
  try(InputStream fi=Files.newInputStream(Path.of(args[0]));InputStream bi=new BufferedInputStream(fi);InputStream gzi=new GzipCompressorInputStream(bi);TarArchiveInputStream archiveStream=new TarArchiveInputStream(gzi)) {
   TarArchiveEntry entry;
   while((entry=archiveStream.getNextEntry())!=null) {
    if(!archiveStream.canReadEntryData(entry))throw new AssertionError("unreadable");
    byte[] bytes=archiveStream.readAllBytes();
    String table=new String(bytes,java.nio.charset.StandardCharsets.UTF_8);
    if(!table.startsWith("#<METADATA>SAMPLE=HCC1395_HCC1395_WES\ncontext\trev_comp\t"))throw new AssertionError("metadata/header");
    if(table.lines().count()!=66)throw new AssertionError("row count");
    System.out.println(entry.getName()+" "+bytes.length+" bytes; GATK archive-reader stream accepted");
   }
  }
 }
}
