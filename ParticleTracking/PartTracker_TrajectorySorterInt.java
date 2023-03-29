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

public class PartTracker_TrajectorySorterInt implements PlugInFilter{
	
	int nFrames;
	int frameSizeX;
	int frameSizeY;
	int nTraj;
	int nPts;
	
	ResultsTable rt;
	
	ImagePlus imp;
	ImageProcessor ip;
	
	Traj[] trajectories;
	ArrayList<int[]>[] ptInT;
	
	boolean canceled;
	String lineSep;
	
	boolean selected[];
	
	// parameters
	
	int minLength;
	int maxLength;
	int nSelecTraj;
	
	public int setup(String arg, ImagePlus imp){
		this.imp = imp;
		ip = imp.getProcessor();
        return DOES_ALL;
	}
	
	public void run(ImageProcessor ip) {
		
		load();
		
		settings();
		
		if(!canceled){
			
			compute();
			
			save();
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
		int ci = rt.getColumnIndex("totInt");
		int cframe = rt.getColumnIndex("frame");
		int ctraj = rt.getColumnIndex("traj");
		int rtlength = rt.getCounter();
		
		if(ctraj==-1){
			Trajectory_Linker tl = new Trajectory_Linker();
			tl.setup("",imp);
			tl.run(ip);
			anal =  new Analyzer(imp); 
			rt = anal.getResultsTable(); 
			cx = rt.getColumnIndex("y");   // Les gars qui ont fait Mosaic, c'est des débiles : ils ont inversé x et y
			cy = rt.getColumnIndex("x"); 
			cs = rt.getColumnIndex("size");
			ci = rt.getColumnIndex("totInt");
			cframe = rt.getColumnIndex("frame");
			ctraj = rt.getColumnIndex("traj");
			rtlength = rt.getCounter();
		}
		
		nTraj = (int) rt.getValueAsDouble(ctraj,rtlength-1);
		nPts = rtlength;
		
		trajectories = new Traj[nTraj];
		int n=1; // traj #
		int ctrL=0;
		nFrames = 0;
		
		for(int i=0; i<rtlength; i++){
			nFrames = ((int)(rt.getValueAsDouble(cframe,i)+1)>nFrames)?(int)(rt.getValueAsDouble(cframe,i)+1):nFrames;
			if((int)rt.getValueAsDouble(ctraj,i)==n){
				ctrL++;
			}
			else{
				trajectories[n-1] = new Traj(ctrL);
				n++;
				ctrL=1;
			}
		}
		trajectories[n-1] = new Traj(ctrL);
		
		ptInT = new ArrayList[nFrames];
		for(int t=0; t<nFrames; t++){
			ptInT[t] = new ArrayList<int[]>();
		}
		
		IJ.log("nb of Traj: "+nTraj);
		IJ.log("nb of Frames: "+nFrames);
		IJ.log("nb of Points: "+nPts);
		
		n=1; // traj #
		ctrL=0; // current index in the trajectory
		
		int[] pointer;
		
		for(int i=0; i<rtlength; i++){
			//IJ.log("i= "+i);
			IJ.showProgress((double)i/rtlength);
			pointer = new int[2];
			pointer[0] = -1+(int)rt.getValueAsDouble(ctraj,i); //traj nb
			pointer[1] = ctrL; // index in traj
			
			trajectories[pointer[0]].x[ctrL] = rt.getValueAsDouble(cx,i);
			trajectories[pointer[0]].y[ctrL] = rt.getValueAsDouble(cy,i);
			trajectories[pointer[0]].t[ctrL] = (int) rt.getValueAsDouble(cframe,i);
			trajectories[pointer[0]].size[ctrL] = (int) rt.getValueAsDouble(cs,i);
			trajectories[pointer[0]].totInt[ctrL] = (int) rt.getValueAsDouble(ci,i);

			ptInT[(int) rt.getValueAsDouble(cframe,i)].add(pointer);
			
			if(i==(rtlength-1)){ctrL=0;}
			else if((int)rt.getValueAsDouble(ctraj,i+1)==n){
				ctrL++;
			}
			else{
				n++;
				ctrL=0;
			}
		}
		
	}
	
	private void settings(){
		
		canceled = false;
		
		GenericDialog gd = new GenericDialog("Postprocessing settings...");
		gd.addNumericField("minimum_trajectory_length (frames)",30,0);
		
		// show the dialog and quit, if the user clicks "cancel"
		gd.showDialog();
		if (gd.wasCanceled()) {
			canceled=true;
			return;
		}
		else {
			minLength = (int) gd.getNextNumber();
			
		}
	}
	
	private void compute(){
		
		Traj traj;
		
		selected = new boolean[trajectories.length];
		
		Roi activeRoi = imp.getRoi();
		if(activeRoi == null){
			activeRoi = new Roi(0,0,frameSizeX,frameSizeY);
		}
		else if(!activeRoi.isArea()){
			activeRoi = new Roi(0,0,frameSizeX,frameSizeY);
		}
		
		maxLength =0;
		nSelecTraj=0;
		
		for(int i=0; i<trajectories.length; i++){
			
			traj = trajectories[i];
			
			selected[i]=(traj.length>minLength) && (activeRoi.contains((int)traj.x[0],(int)traj.y[0]));
			
			nSelecTraj+=(selected[i])?1:0;
			
			if(selected[i] && traj.length>maxLength){
				maxLength = traj.length;
			}
		}
		
	}
	
	private void save(){
		
		String fileName = "ReogTraj_"+imp.getTitle();
		SaveDialog sd = new SaveDialog("saveName",fileName,".txt");
		String dir = sd.getDirectory();
		fileName = sd.getFileName();
		
		Traj traj;
		lineSep = System.getProperty("line.separator");
		String buffer = "t";
		for(int i=0; i<trajectories.length; i++){
			if(selected[i]){
				buffer+="\tmeanInt_Traj_"+i;
			}
		}
		
		try {
			
			FileWriter file = new FileWriter(dir+fileName);
			
			
			buffer += lineSep;
			file.write(buffer);
			
			for (int t=1;t<maxLength+1; t++){
				
				buffer = ""+t;
				
				for(int i=0; i<trajectories.length; i++){
					if(selected[i]){
						traj=trajectories[i];
						if(traj.length>=t){
							buffer+="\t"+(traj.totInt[t-1]/traj.size[t-1]);
						}
						else{
							buffer+="\t";
						}
					}
				}
				
				IJ.showProgress((double)t/maxLength);
				buffer += lineSep;
				file.write(buffer);
				
			}
			file.close();
		} catch (Exception e){
			IJ.log("Erreur doSaveBrown --> "+e.getMessage());
			IJ.log("Erreur doSaveBrown --> "+e.getCause());
			IJ.log("Erreur doSaveBrown --> "+e.getLocalizedMessage());
		} 	
		IJ.showStatus("Done");
		
	}
	
	private int min(int a,int b){
		
		return (a>b)?b:a;
		
	}
	private int max(int a,int b){
		
		return (a<b)?b:a;
		
	}
	
	
}
