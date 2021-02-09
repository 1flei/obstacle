# obstacle
Detecting implicit obstacle from trajectories

This project aims to detect implicit obstacle from trajectory datasets.
It has the following dependencies:

- c++-17
- boost
- mongodb and mongo-cxx-driver
- function_ref (already included)
- hnswlib (already included)

In order to build the project, you should have a c++ compiler with c++-17. We use g++-8 in our experiment, and you can use the following commands to install a g++8 alternative on Ubuntu:

```
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt-get update
sudo apt-get install g++-8
sudo apt-get install gcc-x (x=version)
```

Make sure you set CXX to g++-8 before running the building command by 

```
export CXX=g++-8
```

Before buliding, you should also properly install the boost library and mongodb. After that, you can run

```
mkdir build
cd build
cmake ..
make -j
```

to build the project.

To run the experiment, you can use the following command to run the corresponding experiments. 
```
python ../run_scriptes/run_kde_anomaly_test.py
python ../run_scriptes/run_density_benchmark.py
python ../run_scriptes/run_obst_detection.py
```

Here the python code is just like a scirpts to invoke the corresponding executable files. To change the parameter settings, you can either modify the python code to directly set the parameter to the executable files.
Using the --help parameter to the executable file such as 
```
./obst_test --help
```
can give you the explanation of the parameter for the corresponding program.


To visualize the experiment results, you should first collect the data to the res folder (which is the default output folder). The visualization is done by python and matplotlib. You may need matplotlib to plot the figures.
The plot scripts is at ./plot. Running the corresponding scripts can reproduce the figures used in the paper. 
