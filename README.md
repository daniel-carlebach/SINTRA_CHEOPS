# SINTRA_CHEOPS

This package contains the work that I did in the summer of 2024 as part of the Summer Undergraduate Research in Engineering (SURE) Program at the University of Michigan. My work is part of the Space Debris Identification and Tracking (SINTRA) group. My goal was to analyze images produced by the Characterizing Exoplanets Satellite (CHEOPS) to look for evidence of space debris and collisions between space debris. The progression of my work was first looking for streaks in the images produced by resident space objects (RSOs), then determining which of these streaks were produced by known RSOs, and finally which streaks were likely produced by previously unknown RSOs.

This package contains functions to perform the following 4 tasks:
* Download CHEOPS images and TLE tracking data
* Identify streaks in CHEOPS images
* Determine when RSOs pass within CHEOPS' field of view
* Match streaks in CHEOPS images to transiting RSOs.
