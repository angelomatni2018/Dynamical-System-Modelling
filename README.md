# Dynamical-System-Modelling
The models here are exercises in differential equation solving and analysis, employing several numerical methods of varying complexity, fine tuning them to the specific situations at hand.

## Background to the Mathematics Used

### Numerical Solving Methods
Numerical estimation and smoothing techniques such as the below were used:
* Gradients, Hessians
* Euler
* Runge-Kutta (Order 2 and 4 applied)
* Time-step control, Richardson extrapolation

## Application 1: Bird Flocks

### Objective
A pair of differential equations were analyzed that describe the changes in position and velocity of each bird in a flock, or each entity of a swarm. Numerical methods such as Euler's and Runge-Kutta were employed to understand the behavior of the flock as time went on, given various starting conditions. 

### Implementation
Before behaviour could be analyzed, a reliable method of analysis was found. Euler's approach worked, but had rather large errors. Runge Kutta was far more successful, especially higher order forms of the algorithm. The file ode_solver employs these strategies alongside time-step control and Richardson extrapolation. It then compares them to a MATLAB built in ode solver: ode45.

Once solving is in place, analysis of each bird's trajectory is captured and concisely graphed to show the flock's behaviour. Different seeds are used to vary the initial position and velocity conditions of each bird in the flock.  

### Results
RK4 competes well with ode45 in terms of accuracy and error, so it was used to analyze the behaviour of the flock. It is apparent, as long as initial position and velocity conditions aren't too extreme, that the flock synchronizes and begins flocking in the truest sense of the word. An example of this is shown below:

![](https://github.com/angelomatni2018/Dynamical-System-Modelling/blob/master/bird_flocks/graphs/ang_1_seed17.png?raw=true) | ![](https://github.com/angelomatni2018/Dynamical-System-Modelling/blob/master/bird_flocks/graphs/traj1_seed17.png?raw=true)
