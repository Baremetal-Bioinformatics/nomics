import java.util.*;
public class EmOracle {
 public static void main(String[] args) {
  int ref=1; double[] pseudo=new double[12];Arrays.fill(pseudo,.1);pseudo[1]=pseudo[5]=0;double[] p=pseudo.clone();
  int[][] rows={{30,2,2,0},{30,2,0,0},{100,50,23,0}};
  int it=0;
  while(true) {
   double[] design=new double[12], singles=new double[12], refs=new double[12];
   for(int[] r:rows) {double[] v=Oracle.resp(ref,r[3],r[0],r[1],r[2],p);for(int s=0;s<12;s++)design[s]+=v[s];}
   for(int alt=0;alt<4;alt++)if(alt!=ref)for(int o=0;o<2;o++) {
    double[] effective=new double[12];
    for(int d=1;d<=200;d++) {double[] v=Oracle.resp(ref,alt,d,1,o==0?1:0,p);double count=d==20&&alt==0&&o==0?7:0;for(int s=0;s<12;s++)effective[s]+=count*v[s];}
    for(int s=0;s<12;s++)singles[s]+=effective[s];
   }
   for(int d=1;d<=200;d++) {double[] v=Oracle.resp(ref,ref,d,0,0,p);double count=d==30?1000:0;for(int s=0;s<12;s++)refs[s]+=count*v[s];}
   double[] next=new double[12]; double sum=0;
   for(int s=0;s<12;s++) {next[s]=((design[s]+singles[s])+refs[s])+pseudo[s];sum+=next[s];}
   for(int s=0;s<12;s++)next[s]/=sum;
   double dist=0;for(int s=0;s<12;s++){double v=p[s]-next[s];dist+=v*v;}
   p=next;it++;
   System.out.print("iteration "+it);for(double v:p)System.out.printf(" %016x",Double.doubleToLongBits(v));System.out.println();
   if(Math.sqrt(dist)<=1e-4||it>=20)break;
  }
 }
}
