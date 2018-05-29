# Moving-edge-detection-via-heat-flow-analogy
![Alt text](https://github.com/YuhaoYeSteve/Moving-edge-detection-via-heat-flow-analogy/raw/master/Gif/lab.gif)
Video from: http://www.gait.ecs.soton.ac.uk
181 Frames in total



![Alt text](https://github.com/YuhaoYeSteve/Moving-edge-detection-via-heat-flow-analogy/raw/master/Gif/gait_moving_edge.gif)

Final Result

6 Main Process:

1. Read in three consecutive frames from video
    --Frame(n-1) and frame(n+1) help reference frame(n) find edge  
    --So, processing start from the third frame 

2. Anisotropic Diffusion
    -- AD Vs Gaussian Smoothing(both remove noise)
    -- AD can keep edges(sharpen)

3. Sobel Edge Detection

4. Isotropic heat flow in temporal domain
    -- Heat flow out=Moving edges

5. Non-maximum Suppression
    --Thinning all edges to 1 pixel wide by 
      comparing with local maximum along 
      direction of gradient

6. Hysteresis Thresholding
   -- Filter out false edges & Keep connected “weak” edges
    -- Get rid of background noise  




