package ParticleTracking;

import java.lang.*;
import java.util.*;
import java.io.*;
import java.awt.*;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.frame.*;
import ij.plugin.filter.*;
import ij.plugin.*;
import ij.io.*;
import ij.text.*;
import ij.measure.*;

public class Trajectory_Linker implements PlugInFilter{
	
	int nFrames;
	int frameSizeX;
	int frameSizeY;
	
	ResultsTable rt;
	
	ImagePlus imp;
	ImageProcessor ip;
	
	Point[] pts;
	ArrayList<Integer> PtInT[];
	
	ArrayList<Double> StaticIntensities[];
	ArrayList<Double> DynamicIntensities[];
	
	ArrayList<ArrayList<Point>> trajectories;
	
	boolean canceled;
	String lineSep;
	
	double maxDist;
	
	double thresholdSize; // to distinguish between clamps and single cells
	double[] avgSize; // avg size of a part in a given trajectory
	
	int minLength;
	int maxLength;
	boolean[] clamp;
	boolean[] swimmer;
	double minRad;
	int posFiltLength;
	
	
	boolean velocityFits;
	int velocityFitLength;
	double maxDistVelFit;

	boolean intON;
	
	boolean pairCorr;
	double[] Gr;
	double samplingDist;
	
	boolean printOption;
	
	public int setup(String arg, ImagePlus imp){
		this.imp = imp;
		ip = imp.getProcessor();
        return DOES_ALL;
	}
	
	public void run(ImageProcessor ip) {
		
		load();
		
		settings();
		
		if(!canceled){
			
			trajectoryLinking();
			
			saveTrajectories();
			
			if(printOption){
				plotTrajectories();
			}
			
		}
	}// run
	
	private void load(){
		
		// mise en condition
		
		frameSizeX = imp.getWidth();
		frameSizeY = imp.getHeight();
		lineSep = System.getProperty("line.separator");
		
		// lecture de la table de résultats
		Analyzer anal =  new Analyzer(imp); 
		rt = anal.getResultsTable(); 
		int cx = rt.getColumnIndex("y");   // Les gars qui ont fait Mosaic, c'est des débiles : ils ont inversé x et y
		int cy = rt.getColumnIndex("x"); 
		int cs = rt.getColumnIndex("size"); 
		int cframe = rt.getColumnIndex("frame");
		int cInt = rt.getColumnIndex("totInt");
		//int cscore = rt.getColumnIndex("NPscore");
		int rtlength = rt.getCounter();
		
		nFrames = (int) rt.getValueAsDouble(cframe,rtlength-1)+1;
		
		IJ.log("nb de Frames: "+nFrames);
		
		pts = new Point[rtlength];
		PtInT = new ArrayList[nFrames];
		for(int t=0; t<nFrames; t++){
			PtInT[t] = new ArrayList<Integer>();
		}
		int n=1; // frame #

		intON =  cInt!=ResultsTable.COLUMN_NOT_FOUND;

		if(intON) {
			for(int i=0; i<rtlength; i++) {
				pts[i] = new Point((int) rt.getValueAsDouble(cframe, i), rt.getValueAsDouble(cx, i), rt.getValueAsDouble(cy, i), rt.getValueAsDouble(cs, i), rt.getValueAsDouble(cInt, i));
				PtInT[(int) rt.getValueAsDouble(cframe, i)].add(i);
			}
		}
		else{
			for(int i=0; i<rtlength; i++) {
				pts[i] = new Point((int) rt.getValueAsDouble(cframe, i), rt.getValueAsDouble(cx, i), rt.getValueAsDouble(cy, i), rt.getValueAsDouble(cs, i));
				PtInT[(int) rt.getValueAsDouble(cframe, i)].add(i);
			}
		}
	}
	
	private void settings(){
		
		canceled = false;
		
		GenericDialog gd = new GenericDialog("TrajLinking settings...");
		gd.addNumericField("maximum_Displacement:",10.0,1);
		gd.addCheckbox("print",true);
		
		// show the dialog and quit, if the user clicks "cancel"
		gd.showDialog();
		if (gd.wasCanceled()) {
			canceled=true;
			return;
		}
		else {
			maxDist = (double) gd.getNextNumber();
			printOption = gd.getNextBoolean();
			maxDist*=maxDist;
		}
	}
	
	private void trajectoryLinking(){
		
		trajectories = new ArrayList<ArrayList<Point>>();
		
		double[][] distanceMap;
		int n,m;
		ArrayList<Integer> inT, inTpOne;
		double x0,y0,x1,y1;
		
		double[] minN, minM;
		int[] minDn, minDm;
		
		for(int t=0; t<nFrames-1; t++){
			
			IJ.showProgress((double)t/(nFrames-1));
			
			n = PtInT[t].size();
			m = PtInT[t+1].size();
			
			if(n==0 || m==0){
				continue;
			}
			
			inT = PtInT[t]; // list of the points in t
			inTpOne = PtInT[t+1]; // list of the points in t+1
			
			distanceMap = new double[n][m];
			
			minN=new double[n];
			minM=new double[m];
			minDn=new int[n];
			minDm=new int[m];
			
			for(int j=0; j<m; j++){
				minM[j] = 1e37;
			}
			
			for(int i=0; i<n; i++){
				minN[i] = 1e37;
				
				x0 = pts[inT.get(i)].getX();
				y0 = pts[inT.get(i)].getY();
				
				for(int j=0; j<m; j++){
					x1 = pts[inTpOne.get(j)].getX();
					y1 = pts[inTpOne.get(j)].getY();
					distanceMap[i][j] = (x1-x0)*(x1-x0) +(y1-y0)*(y1-y0);
					
					if(minN[i]>distanceMap[i][j]){
						minN[i]=distanceMap[i][j];
						minDn[i] = j;
					}
					if(minM[j]>distanceMap[i][j]){
						minM[j]=distanceMap[i][j];
						minDm[j] = i;
					}
				}
			}
			
			for(int i=0; i<n; i++){
				if(minDm[minDn[i]]==i && distanceMap[i][minDn[i]]<maxDist){ // i and minDn[i] are reciprocal nearest match and they are not too far away
					
					if(pts[inT.get(i)].getTrajNb()==0){
						ArrayList<Point> newTraj = new ArrayList<Point>();
						pts[inT.get(i)].setTrajNb(trajectories.size()+1);
						pts[inTpOne.get(minDn[i])].setTrajNb(trajectories.size()+1);
						newTraj.add(pts[inT.get(i)]);
						newTraj.add(pts[inTpOne.get(minDn[i])]);
						trajectories.add(newTraj);
					}
					else{
						pts[inTpOne.get(minDn[i])].setTrajNb(pts[inT.get(i)].getTrajNb());
						trajectories.get(pts[inT.get(i)].getTrajNb()-1).add(pts[inTpOne.get(minDn[i])]);
					}
					
				}
				else if(pts[inT.get(i)].getTrajNb()==0){ // taking care of left loners in frame t (should serve only at frame 0)
					ArrayList<Point> newTraj = new ArrayList<Point>();
					pts[inT.get(i)].setTrajNb(trajectories.size()+1);
					newTraj.add(pts[inT.get(i)]);
					trajectories.add(newTraj);
				}
				
			}
			for(int j=0; j<m; j++){
				if(pts[inTpOne.get(j)].getTrajNb()==0){ // this point in t+1 was no match for anybody
					ArrayList<Point> newTraj = new ArrayList<Point>();
					pts[inTpOne.get(j)].setTrajNb(trajectories.size()+1);
					newTraj.add(pts[inTpOne.get(j)]);
					trajectories.add(newTraj);
				}
				
			}
			
		}
		
	}
	
	private void saveTrajectories(){
		
		String fileName = "TrajTracked_"+imp.getTitle();
		SaveDialog sd = new SaveDialog("Fichier Part",fileName,".xls");
		String dir = sd.getDirectory();
		fileName = sd.getFileName();
		
		if(fileName==null){return;}
		
		Point pt;
		ArrayList<Point> traj;

		String buffer ="traj\tframe\tx\ty\tsize"+((intON)?"\ttotInt":"");
		
		try {
			FileWriter file = new FileWriter(dir+fileName);
			
			buffer += lineSep;
			file.write(buffer);
			
			for (int tj=0;tj<trajectories.size(); tj++){
				
				traj = trajectories.get(tj);
				
				for(int i=0; i<traj.size(); i++){
					
					pt = traj.get(i);
					
					buffer = ""+(tj+1)+"\t"+pt.getT()+"\t"+pt.getY()+"\t"+pt.getX()+"\t"+pt.getSize()+((intON)?("\t"+pt.getTotInt()):""); // x and y inverted as in Mosaic old version
					buffer += lineSep;
					
					file.write(buffer);
				}
				
				IJ.showProgress((double)tj/trajectories.size());
				
			}
			file.close();
		} catch (Exception e){
			IJ.log("Erreur doSaveBrown --> "+e.getMessage());
			IJ.log("Erreur doSaveBrown --> "+e.getCause());
			IJ.log("Erreur doSaveBrown --> "+e.getLocalizedMessage());
		} 	
		IJ.showStatus("Done");
		
		IJ.open(dir+fileName);
		
	}
	
	private void plotTrajectories(){
		IJ.run("Remove Overlay", "");
		OvalRoi trace;
		RoiManager rm;
		rm = RoiManager.getInstance();
		if (rm==null) rm = new RoiManager(); 
		else{rm.reset();
		}
		
		ArrayList<Point> curTraj;
		Point curPt;
		
		ArrayList<Double>[] handelInt;
		
		imp.hide();
		
		double regSize;
		
		Color code = new Color(255, 0, 0);
		
		for(int i=0; i<trajectories.size(); i++){
			
			curTraj = trajectories.get(i);
			
			for(int j=0; j<curTraj.size(); j++){
				curPt = curTraj.get(j);
				regSize = Math.sqrt(curPt.getSize()/3.141592);
				trace = new OvalRoi((curPt.getX())-regSize,(curPt.getY())-regSize,2*regSize,2*regSize);
				trace.setStrokeColor(code);
				imp.setSlice(curPt.getT()+1);
				rm.add(imp,trace,1);
				rm.moveRoisToOverlay(imp);
				rm.reset();
			}
		}
		
		Color[] colors = new Color[10];
		colors[0] = new Color(255,  0,  0);
		colors[1] = new Color(0  ,255,  0);
		colors[2] = new Color(0  ,  0,255);
		colors[3] = new Color(255,255,  0);
		colors[4] = new Color(0  ,255,255);
		colors[5] = new Color(255,  0,255);
		colors[6] = new Color(255,255,255);
		colors[7] = new Color(0  ,  0,  0);
		colors[8] = new Color(128,128,128);
		colors[9] = new Color(128,  0,  0);
		
		for(int i=0; i<trajectories.size(); i++){
			
			curTraj = trajectories.get(i);
			
			//*  trajectories
			PolygonRoi trace2;
			
			float[] x = new float[curTraj.size()];
			float[] y = new float[curTraj.size()];
			
			for(int j=0; j<curTraj.size(); j++){
				curPt = curTraj.get(j);
				x[j] = (float) curPt.getX();
				y[j] = (float) curPt.getY();
			}
			
			
			trace2 = new PolygonRoi(x,y,x.length,Roi.POLYLINE);
			trace2.setStrokeColor(colors[i%10]);
			rm.add(imp,trace2,i);
			
			//*/
		}	
		rm.moveRoisToOverlay(imp);
		rm.show();
		
		imp.setSlice(1);
		imp.show();
		
	}
	
	
	
	private int min(int a,int b){
		
		return (a>b)?b:a;
		
	}
	private int max(int a,int b){
		
		return (a<b)?b:a;
		
	}
	
	private double square(double x){
		return x*x;
	}
	
	private double fitLin(double[] xx, int[] tt, int ii, int d){
		
		double un=0;
		double zt=0;
		double z=0;
		double t2=0;
		double t1=0;
		
		int min= max(ii-d/2,0); 
		int max = min(ii+d/2,xx.length-1);
		
		for(int j=min; j<=max; j++){
			
			un += 1;
			zt += tt[j]*xx[j];
			z += xx[j];
			t1 += tt[j];
			t2 += tt[j]*tt[j];
			
		}
		
		return (un*zt - z*t1)/(un*t2 - t1*t1);
		
	}
	
}
