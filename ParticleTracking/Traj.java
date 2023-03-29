package ParticleTracking;

import java.lang.*;
import java.util.*;
import java.io.*;


import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.frame.*;
import ij.plugin.*;
import ij.io.*;
import ij.text.*;
import ij.measure.*;

class Traj {
	
	int length;
	int trajDur;
	int filteredTrajDur;
	double[] x;
	double[] y;
	double[] msd;
	int[] t;
	int[] Nmsd;
	int[] size;
	double[] totInt;
	
	double[] v;  // velocity of the particle as a fct of time
	double[] vn;  // velocity of the particle as a fct of time
	
	double[] vX;  // velocity of the particle as a fct of time
	double[] vY;  // velocity of the particle as a fct of time
	
	
	double[] R2; // trajectory elongation parameter
	
	boolean[] cover;
	
	private ArrayList<Double> runLengths;
	private ArrayList<double[][]> runs;
	private ArrayList<double[][]> tumbles;

	private int tumbleState[]; // 0=undecided, 1=run, 2=tumble
	
	private double varHigh;
	
	// Constructeurs
	Traj(){
		length = 1;
		x = new double[length];
		y = new double[length];
		t = new int[length];
		size = new int[length];
		totInt = new double[length];
		
	}
	Traj(int l){
		length = l;
		x = new double[length];
		y = new double[length];
		t = new int[length];
		size = new int[length];
		totInt = new double[length];
	}
	
	Traj(double[] x, double[] y, int[] t){
		length = x.length;
		if(y.length !=length || t.length!= length){
			throw new IllegalArgumentException();
		}
		else{
			this.x = x;
			this.y = y;
			this.t = t;
		}
		size = new int[length];
		totInt = new double[length];
	}
	
	// routines
	
	void concatenate(Traj traj2){
		
		int newLength = this.length+traj2.length;
		x = Arrays.copyOf(x, newLength);
		y = Arrays.copyOf(y, newLength);
		t = Arrays.copyOf(t, newLength);
		size = Arrays.copyOf(size, newLength);
		totInt = Arrays.copyOf(totInt, newLength);
		
		System.arraycopy(traj2.x, 0, x, this.length, traj2.length);
		System.arraycopy(traj2.y, 0, y, this.length, traj2.length);
		System.arraycopy(traj2.t, 0, t, this.length, traj2.length);
		System.arraycopy(traj2.size, 0, size, this.length, traj2.length);
		System.arraycopy(traj2.totInt, 0, totInt, this.length, traj2.length);
		
		this.length = newLength;
		
	}
	
	
	
	void computeMsd(){
		
		double r;
		
		trajDur = (t[length-1]-t[0]+1);
		
		msd = new double[(t[length-1]-t[0]+1)];
		Nmsd = new int[(t[length-1]-t[0]+1)];
		
		for(int j=0; j<length; j++){
			
			for(int i=j+1; i<length; i++){
				r = (x[i]-x[j]);
				msd[t[i]-t[j]]+= r*r;
				r = (y[i]-y[j]);
				msd[t[i]-t[j]]+= r*r;
				Nmsd[t[i]-t[j]]++;
			}
			
		}
		if(length==1){
			msd[0]=0;
			Nmsd[0]=0;
		}
		for(int j=1; j<trajDur; j++){
			if(Nmsd[j]!=0){
				msd[j]/=Nmsd[j];
			}
		}
		
	}
	
	double gyrationRadius(){
		
		double xm=0;
		double ym=0;
		double x2=0;
		double y2=0;
		
		trajDur = (t[length-1]-t[0]+1);
		
		for(int j=0; j<length; j++){
			
			xm+=x[j];
			ym+=y[j];
			
			x2+=x[j]*x[j];
			y2+=y[j]*y[j];
			
		}
		
		x2/=length;
		y2/=length;
		xm/=length;
		ym/=length;
		
		return (x2-xm*xm+y2-ym*ym)/(trajDur);
	}
	
	double powerFitMsd(double prec, double lnRes){
		
		double un=0;
		double xy=0;
		double y=0;
		double xx=0;
		double x=0;
		
		double p, lnt,lnMsd;
		lnt = -1;
		
		for(int j=1; j<trajDur; j++){
			
			if(Math.log(j)>lnt+lnRes){ // discard points if sampling in log is too fine
				
				lnt = Math.log(j);
				lnMsd = Math.log(msd[j]);
				
				p = prec/min(msd[j],prec) + 1.0/Nmsd[j];
				
				un += p;
				xy += lnt*lnMsd*p;
				x += lnt*p;
				y += lnMsd*p;
				xx += lnt*lnt*p;
				
			}
			
		}
		
		return (un*xy - x*y)/(un*xx - x*x);
		
	}
	
	void computeVelocity(int avgL){
		
		v = new double[length];
		
		double vx,vy;
		
		for(int i=0; i<length; i++){
			
			vx = fitLin(x,t,i,avgL);
			vy = fitLin(y,t,i,avgL);
			
			v[i] = Math.sqrt(vx*vx+vy*vy);
			
		}
		
	}
	
	void computeVelocities(int avgL){
		
		v = new double[length];
		vX = new double[length];
		vY = new double[length];
		
		double vx,vy;
		
		for(int i=0; i<length; i++){
			
			vx = fitLin(x,t,i,avgL);
			vy = fitLin(y,t,i,avgL);
			
			vX[i]=vx;
			vY[i]=vy;
			v[i] = Math.sqrt(vx*vx+vy*vy);
			
		}
		
	}
	
	void computeSwimParam(int avgL){
		
		vn = new double[length];
		
		double vm;
		
		int min; 
		int max;
		int ct;
		
		for(int i=0; i<length; i++){
			min= max(i-avgL/2,0);
			max = min(i+avgL/2,length-1);
			
			vm = 0;
			
			for(int ii=min; ii<max; ii++){
				vm += Math.sqrt((x[ii+1] - x[ii])*(x[ii+1] - x[ii]) + (y[ii+1] - y[ii])*(y[ii+1] - y[ii]));
			}
			vm/=(t[max]-t[min]);
			
			vn[i] = v[i]/vm;
			
		}
		
		// computing 10% higher
		double vnRun = averageHighest(vn, 0.25);

		// normalisation
		for(int i=0; i<vn.length; i++){
			vn[i]/=vnRun;
		}
		
		varHigh = varianceHighest(vn, 0.25);
		
	}
	
	boolean touchBound(int width,int height, int d){
		
		boolean out = false;
		
		for(int i=0; i<length; i++){
			
			out = (out || x[i]<d || x[i]>(width-d) || y[i]<d || y[i]>(height-d)); 
			
		}
		
		return out;
		
	}
	
	boolean cutBrown(int mTD, double tV, int label){
		
		
		/////////  global check on max speed ////////////////
		
		double vr = averageHighest(v, 0.1);
		
		if(vr<tV){
			if(label!=-1){
				IJ.log("speed pb with traj "+label);
			}
			return false;
		}
		
		////////////// check on parts //////////////////////////
		
		int init = 0;
		boolean slow = (v[0]<tV);
		
		cover = new boolean[length];
		Arrays.fill(cover, Boolean.TRUE);
		
		filteredTrajDur = trajDur;
		
		for(int i = 0; i<length; i++){
			
			if((v[i]<tV)!=slow || i==(length-1)){
				if(slow){ // end of slow portion
					if((t[i] - t[init]+1)>mTD){
						for(int j=init; j<=i; j++){
							cover[j]=false;
						}
						filteredTrajDur -= (t[i] - t[init]+1);
					}
					
				}
				else{ // start of slow portion
					init = i;
				}
				slow = (v[i]<tV);
			}
			
			
		}
		
		return true;
		
	}
	
	int computeTumbleNb(double crit, int minTD, double sevTh, int label){
		
		if(filteredTrajDur<trajDur){
			if(label!=-1){
				IJ.log("tumble too long "+label);
			}
		}
		
		int out=0;
		int tl;
		
		ArrayList<Double> xs = new ArrayList<Double>();
		ArrayList<Double> ys = new ArrayList<Double>();
		double[][] coord;
		double[] dummy;
		tumbleState = new int[length];
		
		int td=0;
		
		boolean tumbly = (vn[0]<crit);
		double severity = 0;
		if(cover[0]){
			td++;
			xs.add(x[0]);
			ys.add(y[0]);
			if(tumbly){
				out++;
				severity += (crit - vn[0])*(crit - vn[0]);
			}
		}
		
		runLengths = new ArrayList<Double>();
		
		runs = new ArrayList<double[][]>();
		tumbles = new ArrayList<double[][]>();
		
		//ListIterator<double[][]> runsIt = runs.listIterator();
		//ListIterator<double[][]> tumblesIt = tumbles.listIterator();
		
		for(int i = 1; i<length; i++){

				if(cover[i-1] && !cover[i]){
					if(td>minTD){
						if(tumbly){
							// add a new tumble trajectory
							coord = new double[2][xs.size()];
							coord[0] = toArray(xs);
							coord[1] = toArray(ys);
							tumbles.add(coord);
							for(int j=1; j<=td; j++){
								tumbleState[i-j]=2;
							}
						}
						else{
							// add a new run trajectory
							coord = new double[2][xs.size()];
							coord[0] = toArray(xs);
							coord[1] = toArray(ys);
							runs.add(coord);
							for(int j=1; j<=td; j++){
								tumbleState[i-j]=1;
							}
						}
					}
					xs = new ArrayList<Double>();
					ys = new ArrayList<Double>();
					td=0;
				}
			
			if((vn[i]<crit)!=tumbly){
				tumbly = (vn[i]<crit);
				
				if(cover[i]){
					if(tumbly){
						out++;
						severity=0;
						if(td>minTD){
							runLengths.add((double)td);
							// add a new run trajectory
							coord = new double[2][xs.size()];
							coord[0] = toArray(xs);
							coord[1] = toArray(ys);
							runs.add(coord);
							for(int j=1; j<=td; j++){
								tumbleState[i-j]=1;
							}
						}
					}
					else if(td<=minTD || severity <= sevTh*(varHigh*td*td)/2){
						if(label!=-1 && severity <= sevTh*(varHigh*td*td)/2){
							IJ.log("severity exclusion "+label);
						}
						out--;
						//tumbles.remove(tumbles.size());
					}
					else{
						// add a new tumble trajectory
						coord = new double[2][xs.size()];
						coord[0] = toArray(xs);
						coord[1] = toArray(ys);
						tumbles.add(coord);
						for(int j=1; j<=td; j++){
							tumbleState[i-j]=2;
						}
					}
				}
				xs = new ArrayList<Double>();
				ys = new ArrayList<Double>();
				td=0;
			}
			
			if(cover[i]){
				td++;
				xs.add(x[i]);
				ys.add(y[i]);
				if(tumbly){
					severity += (crit - vn[i])*(crit - vn[i]);
				}
			}
			
		}
		
		if(cover[length-1]){
			if(tumbly){
				coord = new double[2][xs.size()];
				coord[0] = toArray(xs);
				coord[1] = toArray(ys);
				tumbles.add(coord);
				for(int j=1; j<=td; j++){
					tumbleState[length-j]=2;
				}
				xs = new ArrayList<Double>();
				ys = new ArrayList<Double>();
			}
			else {
				coord = new double[2][xs.size()];
				coord[0] = toArray(xs);
				coord[1] = toArray(ys);
				runs.add(coord);
				xs = new ArrayList<Double>();
				ys = new ArrayList<Double>();
				runLengths.add((double)td);
				for(int j=1; j<=td; j++){
					tumbleState[length-j]=1;
				}
			} 
		}
		
		return out;
	}

	int[] getTumbleState(){
		return tumbleState;
	}
	
	ArrayList<Double> getRunLengths(){
		
		return runLengths;
		
	}
	
	ArrayList<double[][]> getRuns(){
		
		return runs;
		
	}
	ArrayList<double[][]> getTumbles(){
		
		return tumbles;
		
	}

	double[] computeTotCorrFct(){

		double[] out = new double[length];

		double[] ut = new double[4 * (length)];;
		double[] vt = new double[4 * (length)];;

		DoubleFFT_1D fFT1D;

		for(int i = 1; i<length; i++){
			ut[2 * i] = vX[i]/v[i]; // dy,dx
			vt[2 * i] = vY[i]/v[i];
		}

		fFT1D = new DoubleFFT_1D(2*(length));
		fFT1D.complexForward(ut);
		fFT1D.complexForward(vt);

		for (int i = 0; i < 2*(length); i++) {
			ut[2 * i] = ut[2 * i] * ut[2 * i] + ut[2 * i + 1] * ut[2 * i + 1];
			ut[2 * i + 1] = 0;
			vt[2 * i] = vt[2 * i] * vt[2 * i] + vt[2 * i + 1] * vt[2 * i + 1];
			vt[2 * i + 1] = 0;
		}
		fFT1D.complexInverse(ut, true);
		fFT1D.complexInverse(vt, true);
		for (int i = 0; i < length; i++) {
			out[i] += (ut[2 * i] + vt[2 * i])/length;
		}
		for(int i=0; i<length; i++){
			out[i]/=out[0];
		}

		return out;
	}

	double[] computeRotCorrFct(){

		int Tmax = 0;
		double th;

		for(double[][] r:runs){
			Tmax = (Tmax<r[0].length-1)?r[0].length-1:Tmax;
		}

		double[] out = new double[Tmax];
		double[] cnt = new double[Tmax];
		double[] ut;
		double[] vt;

		DoubleFFT_1D fFT1D;

		for(double[][] r:runs){

			if(r[0].length>2) {


				ut = new double[4 * (r[0].length - 1)]; // 0 padding
				vt = new double[4 * (r[0].length - 1)];

				for (int i = 0; i < r[0].length - 1; i++) {

					th = Math.sqrt((r[0][i + 1] - r[0][i])*(r[0][i + 1] - r[0][i])+ (r[1][i + 1] - r[1][i])*(r[1][i + 1] - r[1][i]));
					ut[2 * i] = (r[0][i + 1] - r[0][i])/th; // dy,dx
					vt[2 * i] = (r[1][i + 1] - r[1][i])/th;

				}

				fFT1D = new DoubleFFT_1D(2*(r[0].length - 1));
				fFT1D.complexForward(ut);
				fFT1D.complexForward(vt);

				for (int i = 0; i < 2*(r[0].length - 1); i++) {
					ut[2 * i] = ut[2 * i] * ut[2 * i] + ut[2 * i + 1] * ut[2 * i + 1];
					ut[2 * i + 1] = 0;
					vt[2 * i] = vt[2 * i] * vt[2 * i] + vt[2 * i + 1] * vt[2 * i + 1];
					vt[2 * i + 1] = 0;
				}
				fFT1D.complexInverse(ut, true);
				fFT1D.complexInverse(vt, true);
				for (int i = 0; i < r[0].length - 1; i++) {
					out[i] += (ut[2 * i] + vt[2 * i])/(r[0].length - 1);
					cnt[i]++;
				}
			}

		}
		for(int i=0; i<Tmax; i++){
			out[i]/=((cnt[i]==0)?1:cnt[i]);
			out[i]/=out[0];
		}

		return out;
	}
	
	void filter(int l){
		
		double[] xnew = new double[length];
		double[] ynew = new double[length];
		
		int min,max;
		double xm,ym;
		int ct;
		
		for(int i=0; i<length; i++){
			
			min=0;
			while(t[min]<t[i]-l/2){
				min++;
			}
			max=min;
			while(t[max]<t[i]+l/2 && max<length-1){
				max++;
			}
			
			ct=0;
			xm=0;
			ym=0;
			for(int ii=min; ii<=max; ii++){
				xm+=x[ii];
				ym+=y[ii];
				ct++;
			}
			
			xnew[i]=xm/ct;
			ynew[i]=ym/ct;
		}
		
		x = xnew;
		y = ynew;
		
	}
	
	private double fitLin(double[] xx, int[] tt, int ii, int d){
		
		double un=0;
		double zt=0;
		double z=0;
		double t2=0;
		double t1=0;
		
		int min= max(ii-d/2,0); 
		int max = min(ii+d/2,length-1);
		
		for(int j=min; j<=max; j++){
			
			un += 1;
			zt += tt[j]*xx[j];
			z += xx[j];
			t1 += tt[j];
			t2 += tt[j]*tt[j];
			
		}
		
		return (un*zt - z*t1)/(un*t2 - t1*t1);
		
	}
	
	double averageHighest(double[] m, double pct){
		
		int  n = (int)(m.length*pct);
		if(n == 0){n=1;}
		
		double[] sort = new double[n];
		double curVar,tmp;
		
		for(int i=0; i<m.length; i++){
			
			curVar=m[i];
			
			for(int j=0; j<n; j++){
				
				if(curVar>sort[j]){
					tmp = sort[j];
					sort[j] = curVar;
					curVar = tmp;
				}
				
			}
			
		}
		
		double out=0;
		for(int j=0; j<n; j++){
			
			out+= sort[j];
			
		}
		
		
		return out/n;
	}
	
	double varianceHighest(double[] m, double pct){
		
		// beware this computes the variance in Delta x, were xStart belongs to the 25% Highest
		
		int  n = (int)(m.length*pct);
		if(n == 0){n=1;}
		
		double[] sort = new double[n];
		int[] sortIndex = new int[n];
		double curVar,tmp;
		int tmpInd, curInd;
		
		for(int i=0; i<m.length; i++){
			
			curVar=m[i];
			curInd = i;
			
			for(int j=0; j<n; j++){
				
				if(curVar>sort[j]){
					tmp = sort[j];
					sort[j] = curVar;
					curVar = tmp;
					
					tmpInd = sortIndex[j];
					sortIndex[j] = curInd;
					curInd = tmpInd;
				}
				
			}
			
		}
		
		double out=0,mean=0;
		double r;
		for(int j=0; j<n; j++){
			
			if(sortIndex[j]<m.length-1){
				r = m[sortIndex[j]+1] - m[sortIndex[j]];
			}
			else{
				r = m[sortIndex[j]] - m[sortIndex[j]-1];
			}
			
			out+= r*r;
			mean += r;
		}
		mean/=n;
		
		return (out/n - mean*mean);
	}
	
	double[] sortDecrease(double[] m){
		
		int n=m.length;
		
		double[] out = new double[n];
		
		double d;
		int k;
		
		for(int i=0; i<n; i++){
			out[i]=m[i];
		}
		
		for(int i=0; i<n-1; i++){
			d = out[i];
			k = i;
			for(int j=i+1; j<n; j++){
				if(out[j]>d){
					d = out[j];
					k = j;
				}
			}
			out[k] = out[i];
			out[i] = d;
		}
		
		return out;
		
	}
	
	private double min(double a, double b){
		if(a>b) return b;
		else return a;
	}
	private double max(double a, double b){
		if(a>b) return a;
		else return b;
	}
	private int min(int a, int b){
		if(a>b) return b;
		else return a;
	}
	private int max(int a, int b){
		if(a>b) return a;
		else return b;
	}
	
	private double[] toArray(ArrayList<Double> l){
		double[] out = new double[l.size()];
		for(int i=0; i<l.size(); i++){
			out[i] = l.get(i);
		}
		return out;
	}
	
	/*
	Traj(double[] px, double[] py){
		if(px.length!=py.length){
			length = 0;
			x = new double[length];
			y = new double[length];
		}
		else{
			length = px.length;
			x = px;
			y = py;
		}
	}
	Traj(int l, double[] px, double[] py){
		if(px.length==py.length && px.length==l){
			length = l;
			x = px;
			y = py;
		}
		else{
			length = 0;
			x = new double[length];
			y = new double[length];
		}
	}
	*/
	
}