# ParticleTracking
ImageJ Plugins for Particle Tracking

This is a suite of ImageJ plugins for performing 2D particle tracking.

2023/03/29 First batch of uploads
Point.java
Traj.java
MovieNormalizer_.java
Particle_Tracker.java
Particle_Tracker_totInt.java
Trajectory_Linker.java

File descriptions:

Utility files: 
Point.java: Class for storing coordinates of Each Points
Traj.java: Class for storing Points classified as belonging to the same trajectories 

Preprocessing (optional):
MovieNormalizer_.java: Computes the time average of the movie (mean image) and normalizes each frame by dividing by the mean image

Processing: 
Particle_Tracker.java: group above threshold pixels in Points, computes their position by the intensity weighted centroid method and outputs spatiotemporal position and the size (nb of pixels) of each points

Particle_Tracker_totInt.java: 
Same as Particle_Tracker but: 1. the position of each Point is computed as the average position of the pixels (no intensity weighted centroid)  
                              2. the mean intensity of each Point is an output
                              
Trajectory_Linker.java: connects points into a trajectory following a minimal proximity in the next frame algorithm 

Instructions:
Preparations:
- Upload the ParticleTracking folder in the IJ Plugin directory
- Get jtransforms-2.4.jar from croelmiyn/FourierImageAnalysis/packages and upload to Plugin directory

Execution:
1. MovieNormalizer_ : runs without user inputs

2. Particle_Tracker (or Particle_Tracker_totInt):
- Open movie in ImageJ and decide on a threshold if relevant
- Execute: GUI interface pops up
    Parameters:
    
    min_particle_radius: minimum radius a group of pixel needs to have to be considered a trackable point 
    
    separate: if checked, particles with multiple maxima will be split in several particles, using a Gaussian fit of their intensity profiles
    
    Method for threshold: method for separating Points from background
    
    CutOff: _cutOff_ parameter for the method
      
      pixel counts as Point if:
      -  cutoff_above_background: intensity > median(pixels) + _cutOff_ * standardDeviation(pixels), 
      -  threshold: intensity > _cutOff_
      -  percentile: intensity in the _cutOff_ brighter pixel fraction
      -  AbsThresh_above_bgd: intensity > mean(pixels_in_given_frame) + _cutOff_
 
- Output is a Result_Table containing the detected Points in each frame, which needs to be saved (PartTracked_).

3. Trajectory_Linker
- If FIJI: open the PartTracked_ as a ResultTable (File>Import>Results...)
- Run: parameters:
      Maximum Distance: maximum distance a particle is allowed to travel between 2 frames
      print: if checked, the trajectories and cell positions will be displayed on the movie stack
