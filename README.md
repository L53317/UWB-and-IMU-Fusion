# UWB-and-IMU-Fusion
This is a very high-precision integrated navigation software Kalman Filtering method based on IMU and UWB information with MATLAB. In this project, measurements of the inertial measurement unit (IMU) and Ultra-wideband (UWB) are simulated and used to calculate the positions of the robots moving along a given path. The angle and position errors of the robot are also calculated to show the performance.

## Get started
To run this software, please run the `IMUandUWB.m` in MATLAB.

## Results
The results are given as follows.
### The localization trajectory of the robot
<img src="https://github.com/L53317/UWB-and-IMU-Fusion/results/Relative_dislacements_trajectory.jpg" width = 50% height = 50%/>

### Performance

![KF_Fusion_error](results/KF_Fusion_error.jpg)
![Positio_and_error](results/Position_and_Error_in_axis.jpg)
![States_and_error](results/States_and_error.jpg)
![States_and_variance](results/States_and_variance.jpg)

## References and acknowledgments
Gongmin Yan, Precise Strapdown Inertial Navigation System ([PSINS](http://www.psins.org.cn/kydm)) Toolbox for Matlab, Northwestern Polytechnical University.
