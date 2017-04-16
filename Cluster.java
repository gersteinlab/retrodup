/* Cluster discordant reads extracted from parent genes 
 * And hierachically cluster them to call insertion points
 */


import java.lang.*;
import java.util.*;
import java.io.*;

public class Cluster
{
    // Orientations
    private static int oPM = 10; // +-
    private static int oMP =  1; // -+
    private static int oMM =  0; // --
    private static int oPP = 11; // ++


    private int strand = -1,other_strand = -1;
    public  int strand()       { return strand; }
    public  int other_strand() { return other_strand; }

    private int orientation = -1; // 10,01,00,11 => +-,-+,--,++
    public  int orientation() { return orientation; }

    private String chr = "*",other_chr = "*";
    public  String chromosome()       { return chr; }
    public  String other_chromosome() { return other_chr; }

    private int pos =  -1, other_pos = -1;
    public  int position()       { return pos; }
    public  int other_position() { return other_pos; }

    private String read_name = "";
    public  String read_name() { return read_name; }

    private int read_len = 0;
    public  int read_length() { return read_len; }

    private boolean isGood = true;
    public  boolean isGood() { return isGood; }

    public Cluster(String line,int start,int end)
    {
	StringTokenizer toks = new StringTokenizer(line);
	int mid = (start + end)>>1;
	int one3 = start + (end - start)/3;
	int two3 = start + (end - start)*2/3;

	// 1 Read name
	if (!toks.hasMoreTokens()) return;
	read_name = toks.nextToken();

	// 2 Flag
	if (!toks.hasMoreTokens()) return;
	int flag = Integer.parseInt(toks.nextToken());
	if ((flag & 0x0010) > 0) strand       = 0;
	else                     strand       = 1;
	if ((flag & 0x0020) > 0) other_strand = 0;
	else                     other_strand = 1;
	if (strand >  0 && other_strand >  0) orientation = oPP;
	if (strand == 0 && other_strand >  0) orientation = oMP;
	if (strand == 0 && other_strand == 0) orientation = oMM;
	if (strand >  0 && other_strand == 0) orientation = oPM;

	// 3 Chrom name
	if (!toks.hasMoreTokens()) return;
	chr = toks.nextToken();

	// 4 Position
	if (!toks.hasMoreTokens()) return;
	pos = Integer.parseInt(toks.nextToken());
	if (pos < start || pos > end)  isGood = false;

	// 5 Mapping quality
	if (!toks.hasMoreTokens()) return;
	int qual = Integer.parseInt(toks.nextToken());
	if (qual < 15) isGood = false;

	// 6 CIGAR
	if (!toks.hasMoreTokens()) return;
	toks.nextToken();

	// 7 Other chrom name
	if (!toks.hasMoreTokens()) return;
	other_chr = toks.nextToken();
	if (other_chr.equals("*")) isGood = false;

	// 8 Other position
	if (!toks.hasMoreTokens()) return;
	other_pos = Integer.parseInt(toks.nextToken());
	if (other_chr.equals("=") &&
	    other_pos >= start - 1000 &&
	    other_pos <= end   + 1000) isGood = false;

	// 9 Insert size
	if (!toks.hasMoreTokens()) return;
	int span = Integer.parseInt(toks.nextToken());
	if (span < 0) span = -span;
	
	// 10 Read
	if (!toks.hasMoreTokens()) return;
	String read = toks.nextToken();
	read_len = read.length();
	
	if (other_chr.equals("=")) other_chr = chr;
    }

    private Cluster _next = null;
    public  Cluster next()            { return _next; }
    public  Cluster next(Cluster nxt) { return _next = nxt; }

    public boolean clusterable(Cluster c)
    {
	if (!other_chromosome().equals(c.other_chromosome())) return false;
	int ori = orientation,other_ori = c.orientation();
	if ((ori       == oPM || ori       == oMP) &&
	    (other_ori == oPM || other_ori == oMP)) {
	    if (ori != other_ori) {
		if (ori == oPM && other_ori == oMP) {
		    if (other_pos < c.other_position()) return false;
		    if (pos       < c.position())       return false;
		}
		if (ori == oMP && other_ori == oPM) {
		    if (other_pos > c.other_position()) return false;
		    if (pos       > c.position())       return false;
		}
	    }
	    return true;
	}
	if ((ori       == oPP || ori       == oMM) &&
	    (other_ori == oPP || other_ori == oMM)) {
	    if (ori != other_ori) {
		if (ori == oPP && other_ori == oMM) {
		    if (other_pos > c.other_position()) return false;
		    if (pos       < c.position())       return false;
		}
		if (ori == oMM && other_ori == oPP) {
		    if (other_pos < c.other_position()) return false;
		    if (pos       > c.position())       return false;
		}
	    }
	    return true;
	}
	return false;
    }

    public int scoreOn(Cluster clust)
    {
	int ret = -1;
	int n = 0;
	for (Cluster c1 = this;c1 != null;c1 = c1.next())
	    for (Cluster c2 = clust;c2 != null;c2 = c2.next()) {
		if (!c1.clusterable(c2)) return -1;
		double delta = c2.other_position() - c1.other_position();
		if (delta < 0) delta = -delta;
		if (ret < 0) ret = 0;
		ret += delta;
		n++;
	    }
	if (n > 0) ret /= n;
	return ret;
    }

    public int count() { return count(-1); }
    public int count(int ori)
    {
	Cluster cl = this;
	int ret = 0;
	while (cl != null) {
	    if (ori < 0) ret++;
	    else if (cl.orientation() == ori) ret++;
	    cl = cl.next();
	}
	return ret;
    }
    
    public String toString()
    {
	StringBuffer ret = new StringBuffer();
	ret.append(chromosome()); ret.append(" ");
	ret.append(position());   ret.append(" ");
	ret.append(strand());     ret.append(" ");
	ret.append(other_chromosome()); ret.append(" ");
	ret.append(other_position());   ret.append(" ");
	ret.append(other_strand());     ret.append(" ");
	return ret.toString();
    }

    public void printBed()
    {
	System.out.println("track name=\"Reads\" useScore=1");
	int n = 0, step = 1,score = 0;
	for (Cluster c = this;c != null;c = c.next()) n++;
	if (n < 1000) step = 1000/n;
	for (Cluster c = this;c != null;c = c.next()) {
	    System.out.print("chr" + c.chromosome() + "\t" +
			     c.position() + "\t" +
			     (c.position() + c.read_length()) + "\t" +
			     c.read_name() + "\t" + score + "\t");
	    if (c.strand() == 0) System.out.println("-");
	    else                 System.out.println("+");
	    score += step;
	}
	score = 0;
	for (Cluster c = this;c != null;c = c.next()) {
	    System.out.print("chr" + c.other_chromosome() + "\t" +
			     c.other_position() + "\t" +
			     (c.other_position() + c.read_length()) + "\t" +
			     c.read_name() + "\t" + score + "\t");
	    if (c.other_strand() == 0) System.out.println("-");
	    else                       System.out.println("+");
	    score += step;
	}
    }

    public static void main(String[] args)
    {
	int start = 0,end = 0;
	try {
	    start = Integer.parseInt(args[0]);
	    end   = Integer.parseInt(args[1]);
	} catch (Exception e) {
	    System.err.println("Parsing input arguments failed.");
	    return;
	}

	ArrayList<Cluster> clusters = new ArrayList<Cluster>(100);
	try {
	    BufferedReader stdin =
		new BufferedReader(new InputStreamReader(System.in));
	    String line = "";
	    while ((line = stdin.readLine()) != null) {
		Cluster cl = new Cluster(line,start,end);
		if (cl.isGood()) {
		    clusters.add(cl);
		    //System.out.println(line);
		}
	    }
        } catch (Exception e) {
            System.err.println(e.toString());
        }

	System.out.println("Read " + clusters.size() + " RPs");
	if (clusters.size() > 20000) return;
	

	boolean updated = true;
	while (updated) {
	    int len = clusters.size();
	    if (len%100 == 0) System.out.println("Left " + len + " RPs ...");
	    int best_score = -1, best_i = -1,best_j = -1;
	    for (int i = 0;i < len;i++) {
		Cluster cl1 = clusters.get(i);
		for (int j = i + 1;j < len;j++) {
		    Cluster cl2 = clusters.get(j);
		    int score = cl1.scoreOn(cl2);
		    if (score < 0) continue;
		    if (best_score == -1 || score < best_score) {
			best_score = score;
			best_i    = i;
			best_j    = j;
		    }
		}
	    }
	    if (best_i < 0 || best_j < 0) break;
	    Cluster cl1 = clusters.get(best_i);
	    Cluster cl2 = clusters.get(best_j);
	    if (best_score < 500) {
		while (cl1.next() != null) cl1 = cl1.next();
		cl1.next(cl2);
		clusters.remove(best_j);
		updated = true;
	    } else updated = false;
	}

	int len = clusters.size();
	for (int i = 0;i < len;i++) {
	    Cluster cl = clusters.get(i);
	    int c10 = cl.count(10),c1 = cl.count(1),c11 = cl.count(11),c0 = cl.count(0);
	    int min = 1;
	    int min2 = 10000000;
	    boolean printCluster = false;
	    int s = -1,e = -1;
	    if ((c10 > min && c1 > min) || (c10 > min2 || c1 > min2)) {
		printCluster = true;
		for (Cluster c = cl;c != null;c = c.next()) {
		    int other_pos = c.other_position();
		    if (c.orientation() == oPM &&
			(e < 0 || other_pos < e)) e = other_pos;
		    if (c.orientation() == oMP &&
			(s < 0 || other_pos > s)) s = other_pos;
		}
		System.out.println("Cluster " + cl.count() +
				   " same chr" + cl.other_chromosome() +
				   ":" + s + "-" + e + " " + (e - s));
	    }
	    if ((c11 > min && c0 > min) || (c11 > min2 || c0 > min2)) {
		printCluster = true;
		for (Cluster c = cl;c != null;c = c.next()) {
		    int other_pos = c.other_position();
		    if (c.orientation() == oMM &&
			(e < 0 || other_pos < e)) e = other_pos;
		    if (c.orientation() == oPP &&
			(s < 0 || other_pos > s)) s = other_pos;
		}
		System.out.println("Cluster " + cl.count() +
				   " different chr" +
				   cl.other_chromosome() + ":" + s +
				   "-" + e + " " + (e - s));
	    }
	    if (printCluster) {
		for (Cluster c = cl;c != null;c = c.next())
		    System.out.println(c.toString());
		//cl.printBed();
	    }
	}
    }
    
}
