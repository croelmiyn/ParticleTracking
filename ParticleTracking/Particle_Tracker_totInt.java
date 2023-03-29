package ParticleTracking;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import ij.measure.ResultsTable;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.io.FileWriter;
import java.util.ArrayList;

public class Particle_Tracker_totInt implements PlugInFilter{
	
	int nFrames;
	int frameSizeX;
	int frameSizeY;
	
	ResultsTable rt;
	
	ImagePlus imp;
	ImageProcessor ip;
	
	boolean canceled;
	
	
	ArrayList<Point> PtInT[];
	
	// static variables
	double sqrt2pi = Math.sqrt(2*Math.PI);
	String lineSep = System.lineSeparator();
	
	String[] methods = new String[4];
	
	
	// parameters
	double factor;
	double threshSize;
	boolean separateTouching;
	String methodThresh;
	
	
	public int setup(String arg, ImagePlus imp){
		this.imp = imp;
		ip = imp.getProcessor();
		
		nFrames = imp.getStackSize();
		frameSizeX = imp.getWidth();
		frameSizeY = imp.getHeight();
		
        return DOES_ALL;
	}
	
	public void run(ImageProcessor ip) {
		
		settings();
		
		if(!canceled){
			
			compute();
			
			print();
			
		}
	}// run
	
	
	private void settings(){
		
		canceled = false;
		
		threshSize = 8;
		factor = 4.5;
		
		methods[0] = "cutoff_above_background";
		methods[1] = "threshold";
		methods[2] = "percentile";
		methods[3] = "AbsThresh_above_bgd";
		
		GenericDialog gd = new GenericDialog("Tracking settings...");
		gd.addNumericField("min_particle_radius (px):",threshSize,1);
		gd.addCheckbox("separate touching particles",true);
		gd.addChoice("Method for threshold", methods, methods[0]);
		gd.addNumericField("cutOff_for_particle_determination:",factor,1);
		
		// show the dialog and quit, if the user clicks "cancel"
		gd.showDialog();
		if (gd.wasCanceled()) {
			canceled=true;
			return;
		}
		else {
			threshSize = (double) gd.getNextNumber();
			separateTouching = gd.getNextBoolean();
			methodThresh = gd.getNextChoice();
			factor = (double) gd.getNextNumber();
			
		}
	}
	
	private void compute(){
		
		float[][] intensity;
		
		int[] hist;
		int index;
		int ct;
		
		PtInT = new ArrayList[nFrames];
		
		// for lin fit
		double sIC, sC, sII, sI, sOne;
		
		double slope, threshold;
		
		// for particle segmentation
		int[][] partNb;
		int nPart;
		ArrayList<ArrayList<int[]>> partCoords;
		int[] coord = new int[2];
		int[] coord2 = new int[2];
		ArrayList<int[]> hdlPt;
		int[] nghX = new int[4];
		int[] nghY = new int[4];
		
		float[] intList;
		int[] xList;
		int[] yList;
		int[] partList;
		int[] rank;
		double dist;
		
		IJ.log("compiled");
		
		for(int t=0; t<nFrames; t++){
			
			PtInT[t] = new ArrayList<Point>();
			
			imp.setSlice(t+1);
			intensity = ip.getFloatArray();
			
			if(methodThresh == methods[0]){
				
				// make image histogram
				hist = new int[65536];
				
				for(int x=0; x<frameSizeX; x++){
					for(int y=0; y<frameSizeY; y++){
						hist[(int)intensity[x][y]]++;
					}
				}
				
				// compute typical value
				ct=0;
				index=0;
				for(int i=0; i<hist.length; i++){
					if(ct<hist[i]){
						ct = hist[i];
						index=i;
					}
				}
				
				// compute variance around typical value
				sIC=0;
				sC=0;
				sII=0;
				sI=0;
				sOne=0;
				ct=0;
				for(int i=index-10; i<=index+10; i++){ // set 10 as an option ?
					ct+=hist[i];
					sIC+=ct*i;
					sC+=ct;
					sII+=i*i;
					sI+=i;
					sOne+=1;
				}
				slope = (sOne*sIC-sI*sC)/(sII*sOne-sI*sI);
				
				// set the threshold pixel value
				threshold = index + factor * frameSizeX*frameSizeY / (slope*sqrt2pi);
			}
			else if(methodThresh == methods[1]){
				threshold = factor;
			}
			else if(methodThresh == methods[2]){
				// make image histogram
				hist = new int[65536];
				
				for(int x=0; x<frameSizeX; x++){
					for(int y=0; y<frameSizeY; y++){
						hist[(int)intensity[x][y]]++;
						
					}
				}
				
				// compute typical value
				double half = frameSizeX*frameSizeY*(1-factor);
				ct=0;
				index=0;
				while(ct<half && index<65535){
					ct += hist[index];
					index++;
				}
				threshold = index;
			}
			else if(methodThresh == methods[3]){
				// make image histogram
				hist = new int[65536];
				
				for(int x=0; x<frameSizeX; x++){
					for(int y=0; y<frameSizeY; y++){
						hist[(int)intensity[x][y]]++;
					}
				}
				
				// compute typical value
				ct=0;
				index=0;
				for(int i=0; i<hist.length; i++){
					if(ct<hist[i]){
						ct = hist[i];
						index=i;
					}
				}
				
				threshold = index + factor;
			}
			else{
				threshold = factor;
			}
			
			// find all pixels above threshold and group as particles
			partNb = new int[frameSizeX][frameSizeY];
			partCoords = new ArrayList<ArrayList<int[]>>();
			nPart=0;
			
			for(int x=1; x<frameSizeX-1; x++){
				for(int y=1; y<frameSizeY-1; y++){
					
					if(intensity[x][y]>threshold){ // belongs to a particle
						
						coord = new int[2];
						coord[0]=x;
						coord[1]=y;
						
						if(partNb[x-1][y-1]==0 && partNb[x-1][y]==0 && partNb[x-1][y+1]==0 && partNb[x][y-1]==0){ // new part
							
							nPart++;
							partNb[x][y]=nPart;
							hdlPt = new ArrayList<int[]>();
							hdlPt.add(coord);
							partCoords.add(hdlPt);
							
						}
						else{ // new px in old point
							nghX = new int[4];
							nghY = new int[4];
							nghX[0]=x-1;
							nghX[1]=x-1;
							nghX[2]=x-1;
							nghX[3]=x;
							nghY[0]=y-1;
							nghY[1]=y;
							nghY[2]=y+1;
							nghY[3]=y-1;
							
							
							for(int i=0; i<4; i++){
								if(partNb[nghX[i]][nghY[i]]!=0){ // point in a cluster
									if(partNb[x][y]==0){ // add a new pixel in the cluster
										partNb[x][y]=partNb[nghX[i]][nghY[i]];
										partCoords.get(partNb[x][y]-1).add(coord);
										//IJ.log("add part: "+coord[0]+"  "+coord[1]+"  from part: "+nghX[i]+"  "+nghY[i]);
									}
									else if(partNb[x][y]!=partNb[nghX[i]][nghY[i]]){ // cluster merging
										int mergeCl = min(partNb[x][y],partNb[nghX[i]][nghY[i]])-1; // list index
										int destrCl = max(partNb[x][y],partNb[nghX[i]][nghY[i]])-1;
										
										//IJ.log(""+x+"  "+y+"  "+nghX[i]+"  "+nghY[i]+"  "+mergeCl+"    "+destrCl+"  "+nPart);
										//IJ.log(""+x+"  "+y+"  "+nghX[i]+"  "+nghY[i]+"  "+partNb[x][y]+"    "+partNb[nghX[i]][nghY[i]]+"  "+nPart);
										
										hdlPt = partCoords.get(destrCl);
										//IJ.log(""+hdlPt.size());
										for(int k = 0; k<hdlPt.size(); k++){
											coord2 = hdlPt.get(k);
											partNb[coord2[0]][coord2[1]]=mergeCl+1;
											//IJ.log(""+coord2[0]+"  "+coord2[1]);
										}
										
										for(int j=destrCl+1; j<nPart; j++){ // rename other clusters
											hdlPt = partCoords.get(j);
											for(int k = 0; k<hdlPt.size(); k++){
												coord2 = hdlPt.get(k);
												partNb[coord2[0]][coord2[1]]--;
											}
										}
										//IJ.log(""+partCoords.get(mergeCl).size());
										partCoords.get(mergeCl).addAll(partCoords.get(destrCl));
										//IJ.log(""+partCoords.get(mergeCl).size());
										
										
										partCoords.remove(destrCl);
										nPart--;
										
									}
								}
							}// neighbors
							
						}
					} // belongs to a part
					
					
				}
			}
			
			IJ.log(""+threshold+"    "+nPart);
			
			// partCoord is an array of groups of connected px above threshold
			// find groups large enough to be particles, compute center, save coordinates
			double[] x0,y0,norm;
			int ptNb;
			ArrayList<Integer> maxima;
			ArrayList<Integer> subSizes;
			
			int[] partSize;
			
			for(int i=0; i<partCoords.size(); i++){
				
				hdlPt = partCoords.get(i); // convenience of writing for the current gp of px
				
				if(hdlPt.size() > 3.14*threshSize*threshSize){ // is a particle
					
					//IJ.log("--------");
					
					///////// separating if there are two or more points close together
					maxima = new ArrayList<Integer>();
					subSizes = new ArrayList<Integer>(); 
					
					intList = new float[hdlPt.size()];
					xList = new int[hdlPt.size()];
					yList = new int[hdlPt.size()];
					partList = new int[hdlPt.size()];
					
					for(int j=0; j<hdlPt.size(); j++){
						coord = hdlPt.get(j);
						xList[j] = coord[0];
						yList[j] = coord[1];
						intList[j] = intensity[coord[0]][coord[1]];
					}
					
					
					if(separateTouching){
						
						// sort intensities
						float d;
						int curL,dumRank;
						rank = new int[hdlPt.size()];
						for(int j=0; j<xList.length; j++){
							rank[j]=j;
						}
						// I want rank[j] is the index of the jth intensity
						for(int j=0; j<rank.length; j++){
							d = intList[rank[j]];
							curL = j;
							for(int l=j+1; l<rank.length; l++){
								if(intList[rank[l]]>d){
									d = intList[rank[l]];
									curL = l;
								}
							} // rank[curL] is the index of the next highest intensity
							dumRank = rank[j];
							rank[j] = rank[curL];
							rank[curL]=dumRank;
						}
						// I am supposed to have rank[j] as the index of the jth part in intensity order
						
						// find isolated maxima
						ptNb = 0;
						int[] minK;
						double[] minD;
						int nbgPt;
						
						for(int j=0; j<xList.length; j++){
							
							//IJ.log(""+intList[rank[j]]);
							minD = new double[ptNb+1];
							minK = new int[ptNb+1];
							for(int c=1; c<ptNb+1; c++){
								minK[c] = -1;
								minD[c] = xList.length*xList.length;
							}
							for(int k=0; k<j; k++){
								dist = xList[rank[j]] - xList[rank[k]];
								dist*=dist;
								dist+= (yList[rank[j]] - yList[rank[k]])*(yList[rank[j]] - yList[rank[k]]);
								if(dist<3 && dist<minD[partList[rank[k]]]){ //here it is 8 connected neighbors for maxima isolation. before was threshSize*threshSize
									minD[partList[rank[k]]] = dist;
									minK[partList[rank[k]]] = rank[k];
								}
							}
							
							nbgPt=0;
							for(int c=1; c<=ptNb; c++){
								if(minK[c]>=0){
									nbgPt++;
								}
							}
							
							if(nbgPt == 0){ // no neighbors
								ptNb ++;
								partList[rank[j]] = ptNb;
								maxima.add(rank[j]);
								subSizes.add(1);
							}
							else if(nbgPt == 1){ // just one neighbor
								for(int c=1; c<=ptNb; c++){
									if(minK[c]>=0){
										partList[rank[j]] = c;
										subSizes.set(c-1,subSizes.get(c-1)+1);
									}
								}
							}
							else{ // consider merging 
								
								ArrayList<Integer> toBeMerged = new ArrayList<Integer>();
								int largestCluster, curSize, absMin;
								largestCluster = 1;
								
								while(minK[largestCluster]<0){
									largestCluster++;
								}
								//IJ.log("lc "+largestCluster);
								double absMinD = xList.length*xList.length;
								absMin = -1;
								
								for(int c=1; c<=ptNb; c++){
									if(minK[c]>=0){
										
										curSize = subSizes.get(c-1);
										//IJ.log(""+c+" "+curSize);
										if(subSizes.get(largestCluster-1)<curSize){
											largestCluster = c;
										}
										if(curSize < (3.14*threshSize*threshSize)){
											toBeMerged.add(c);
										}
										if(minD[c]<absMinD){
											absMin = minK[c];
											absMinD = minD[c];
										}
										
									}
								}
								for(int l=toBeMerged.size()-1; l>=0; l--){ // merge with the largest around, the cluster nbs are ranked increasing in toBeMerged
									
									int mergedClusterNb = toBeMerged.get(l);
									
									if(mergedClusterNb!=largestCluster){
										subSizes.set(largestCluster-1,subSizes.get(largestCluster-1)+subSizes.get(mergedClusterNb-1));
										subSizes.remove(mergedClusterNb-1);
										maxima.remove(mergedClusterNb-1);
										ptNb--;
										if(mergedClusterNb<largestCluster){
											largestCluster--;
										}
										
										for(int k=0; k<j; k++){
											if(partList[rank[k]]==mergedClusterNb){
												partList[rank[k]] = largestCluster;
											}
											else if(partList[rank[k]]>mergedClusterNb){
												partList[rank[k]]--;
											}
										}
									}
									
								}
								//IJ.log("--");
								
								partList[rank[j]] = partList[absMin]; // point attributed to the nearest merged cluster;
								subSizes.set(partList[rank[j]]-1,subSizes.get(partList[rank[j]]-1)+1);
								
							}
							
							
						} // find isolated maxima
						
						// size check
						partSize = new int[ptNb+1];
						for(int j=0; j<xList.length; j++){
							partSize[partList[j]]++;
						}
						for(int k=1; k<=ptNb; k++){
							//IJ.log(""+partSize[k]+"  "+subSizes.get(k-1));
							if(partSize[k]<3.14*threshSize*threshSize){
								for(int j=0; j<xList.length; j++){
									if(partList[j]==k){partList[j]=0;}
								}
							}
						}
						
					}//if separate touching
					
					else{
						ptNb=1;
						for(int j=0; j<xList.length; j++){
							partList[j] = 1;
						}
						partSize = new int[ptNb+1];
						partSize[1] = xList.length;
					}
					
					
					//// not using the centroid method for each subpoint: modif 2021 09 20
					x0 = new double[ptNb+1];
					y0 = new double[ptNb+1];
					norm = new double[ptNb+1];

					for(int j=0; j<partList.length; j++){
						x0[partList[j]] += (double) (xList[j]);
						y0[partList[j]] += (double) (yList[j]);
						norm[partList[j]] += (double) intList[j];
						//IJ.log(""+partList[i]);
					}
					
					for(int j=1; j<=ptNb; j++){
						if(norm[j]!=0){
							x0[j] /= partSize[j];
							y0[j] /= partSize[j];
							PtInT[t].add(new Point(t,x0[j],y0[j],partSize[j],norm[j]));
							//IJ.log(""+x0[j]+"  "+y0[j]);
						}
					}
				} // is part
				
			} // for putative coords
			
		}//for t
		
		
	} // compute
	
	private void print(){
		
		String fileName = "PartTracked_"+imp.getTitle();
		SaveDialog sd = new SaveDialog("Fichier Part",fileName,".xls");
		String dir = sd.getDirectory();
		fileName = sd.getFileName();
		
		if(fileName==null){return;}
		
		Point pt;
		
		String buffer ="frame\tx\ty\tsize\ttotInt";
		
		try {
			FileWriter file = new FileWriter(dir+fileName);
			
			buffer += lineSep;
			file.write(buffer);
			
			for (int t=0;t<nFrames; t++){
				
				for(int i=0; i<PtInT[t].size(); i++){
					
					pt = PtInT[t].get(i);
					
					buffer = ""+pt.getT()+"\t"+pt.getY()+"\t"+pt.getX()+"\t"+pt.getSize()+"\t"+pt.getTotInt(); // x and y inverted as in Mosaic old version
					buffer += lineSep;
					
					file.write(buffer);
				}
				
				IJ.showProgress((double)t/nFrames);
				
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
	
	private int min(int a, int b){
		return (a<b)?a:b;
		
	}
	
	private int max(int a, int b){
		return (a>b)?a:b;
		
	}
}