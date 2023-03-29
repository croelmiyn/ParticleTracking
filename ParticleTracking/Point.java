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

class Point{

	int t;
	double x,y;
	
	int trajNb;
	
	double size;
	double totInt;
	
	Point(int t, double x, double y){
		this.t = t;
		this.x = x;
		this.y = y;
		trajNb=0;
	}
	
	Point(int t, double x, double y, double s){
		this.t = t;
		this.x = x;
		this.y = y;
		size = s;
		trajNb=0;
	}

	Point(int t, double x, double y, double s, double tot){
		this.t = t;
		this.x = x;
		this.y = y;
		size = s;
		totInt = tot;
		trajNb=0;
	}
	
	double getX(){
		return x;
	}
	
	double getY(){
		return y;
	}
	
	int getT(){
		return t;
	}
	
	int getTrajNb(){
		return trajNb;
	}
	
	void setTrajNb(int n){
		trajNb = n;
	}
	
	void setSize(double s){
		size = s;
	}
	
	double getSize(){
		return size;
	}

	void setTotInt(double s){
		totInt = s;
	}
	double getTotInt(){ return totInt;}
}
	
