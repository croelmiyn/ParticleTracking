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

public class MovieNormalizer_ implements PlugInFilter{
	
	ImagePlus imp;
	ImageProcessor ip;
	
	int nFrames;
	int frameSizeX;
	int frameSizeY;
	
	public int setup(String arg, ImagePlus imp){
		this.imp = imp;
		
		nFrames = imp.getStackSize();
		frameSizeX = imp.getWidth();
		frameSizeY = imp.getHeight();
		
        if (arg.equals("about"))
            {return DONE;}
        return DOES_ALL;
	}
	
	public void run(ImageProcessor ip) {
		
		float[][] avgI = new float[frameSizeX][frameSizeY];
		
		imp.hide();
		
		// average
		for(int i = 0; i<nFrames; i++){
			
			IJ.showProgress((double)i/nFrames);
			imp.setSlice(i+1);
			add(avgI,ip.getFloatArray());
			
		}
		
		for(int x=0; x<frameSizeX; x++){
			for(int y=0; y<frameSizeY; y++){
				avgI[x][y]/=nFrames;
				
			}
		}
		
		
		// min and max
		float[][] tempI = new float[frameSizeX][frameSizeY];
		float min = 1; 
		float max = 0;
		float tp;
		
		for(int i = 0; i<nFrames; i++){
			IJ.showProgress((double)i/nFrames);
			
			imp.setSlice(i+1);
			
			tempI = ip.getFloatArray();
			
			for(int x=0; x<frameSizeX; x++){
				for(int y=0; y<frameSizeY; y++){
					tp = tempI[x][y]/avgI[x][y];
					max = (max<tp)?tp:max;
					min = (min>tp)?tp:min;
				}
			}
		}
		
		// 
		ImagePlus imp2 = IJ.createImage("Result_"+imp.getTitle(),"8-bits", frameSizeX, frameSizeY, nFrames);
		ImageProcessor ip2 = imp2.getProcessor();
		
		max-=min; // speed up
		
		for(int i = 0; i<nFrames; i++){
			IJ.showProgress((double)i/nFrames);
			
			imp.setSlice(i+1);
			imp2.setSlice(i+1);
			tempI = ip.getFloatArray();
			
			for(int x=0; x<frameSizeX; x++){
				for(int y=0; y<frameSizeY; y++){
					tp = tempI[x][y]/avgI[x][y];
					ip2.putPixel(x,y,(int)(255*(tp-min)/max));
				}
			}
		}
		
		imp.show();
		imp2.show();
	}
	
	private void add(float[][] a, float[][] b){
		
		for(int x=0; x<a.length; x++){
			for(int y=0; y<a[0].length; y++){
				a[x][y]+=b[x][y];
			}
			
		}
		
	}
}