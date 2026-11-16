import java.util.*;
import java.util.stream.*;
public class HashMapOracle{
 public static void main(String[] args){var keys=new ArrayList<String>();for(int i=0;i<32;i++){var key=new StringBuilder();for(int bit=0;bit<5;bit++)key.append((i&(1<<bit))==0?"Aa":"BB");keys.add(key.toString());}for(int i=0;i<100;i++)keys.add("contig"+i);keys.addAll(Arrays.asList("AaAaAaAaAaAa","BB",keys.get(0),keys.get(20)));var groups=IntStream.range(0,keys.size()).boxed().collect(Collectors.groupingBy(keys::get));for(var group:groups.values())System.out.println(group.get(0));}
}
