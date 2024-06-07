# Compilation guide
## MPBoot
### Downloading source code
You can clone the source code from GitHub with:

`git clone https://github.com/diepthihoang/mpboot.git`
### Compiling under Linux
1. Create folder **build** outside folder **mpboot**.
2. Open Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../mpboot -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.  
  
The compiler will generate an executable file named **mpboot-avx** 
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy**, run command:  
`./mpboot-avx -s example.phy`

### Compiling under Mac OS X
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../mpboot -DIQTREE_FLAGS=avx`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.  
  
The compiler will generate an executable file named **mpboot-avx** 
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy**, run command:  
`./mpboot-avx -s example.phy`

### Compiling under Windows
* Requirements:  
  * cmake version >= 3.21
  * TDM-GCC
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake -G "MinGW Makefiles" -DIQTREE_FLAGS=avx ../mpboot`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture.  
> Due to having conflicts with **Vectorization**, please not using **Clang** to configure source code.

5. Run command `mingw32-make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate two executable files named **mpboot-avx** and **mpboot-avx-click**
> In case of running MPBoot on SSE architecture, the executable files are named **mpboot** and **mpboot-click**.

6. To analyst file **example.phy**:
* Run command `./mpboot-avx -s example.phy`
* Or:
    1. Open **mpboot-click** by double click on it
    2. Press **y** to start enter command
    3. Type **-s** and press **Enter**
    4. Press **e** to continue entering command
    5. Type **example.phy** and press **Enter**
    6. Press **y** to finish command


  
## MPBoot-MPI
### Downloading source code
You can clone the source code from GitHub with:

`git clone https://github.com/diepthihoang/mpboot.git`  
Switch to branch which contains MPBoot-MPI source code:
* `git checkout mpboot-mpi-sync` for synchronous version
* `git checkout mpboot-mpi-async` for asynchronous version

### Compiling under Linux
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../source -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx**
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy** with 2 processes, run command:  
`mpirun -np 2 ./mpboot-avx -s example.phy`
> Option **np** specifies the number of processes used to run MPBoot-MPI

### Compiling under Mac OS X
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../mpboot -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes. Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx**
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy** with 2 processes, run command:  
`mpirun -np 2 ./mpboot-avx -s example.phy` 
> Option **np** specifies the number of processes used to run MPBoot-MPI

### Compiling under Windows
* Requirements:  
  * cmake version >= 3.21
  * TDM-GCC
  * MSMPI
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake -G "MinGW Makefiles" -DIQTREE_FLAGS=mpiavx ../mpboot`
> Replace mpiavx by mpisse4 in above command if you decide to run MPBoot on SSE architecture.  
> Due to having conflicts with **Vectorization**, please not using **Clang** to configure source code.

5. Run command `mingw32-make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx** 
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.
  

6. To analyst file **example.phy** with 2 processes, run command:  
`mpiexec -n 2 ./mpboot-avx -s example.phy`
> Option **n** specifies the number of processes used to run MPBoot-MPI



## MPBoot-SA
MPBoot-SA uses simulated annealing to improve tree-searching phase 
### Downloading source code and building
You can clone the source code and switch to MPBoot-SA version with:

`git clone https://github.com/a2dalek/mpboot.git`

`git checkout sa-reheat-new-global`

### Supported options
To run MPBoot-SA, please add option **-pure_sa** to command. For example, the command used to analyst file **example.phy** by MPBoot-SA with all default options is:

`mpboot-avx -s example.phy -pure_sa`

Currently supported options:
  * `-cooling_schedule <integer>`: Used to choose cooling schedule (defaulted to 2). There are 4 installed schedules, numbered from 0 to 3:
    * EXPONENTIAL MULTIPLICATIVE COOLING
    * EXPONENTIAL ADDITIVE COOLING
    * LINEAR ADDITIVE COOLING
    * LINEAR MULTIPLICATIVE COOLING
  * `-start_temp <float>`: Used to specify starting temperature
  * `-final_temp <float>`: Used to specify the minimum temperature

If not be specified, the default values of **start_temp** and **final_temp** will be determined base on **cooling_schedule**:
  * EXPONENTIAL MULTIPLICATIVE COOLING: start_temp = 0.005, final_temp = 0.00025
  * EXPONENTIAL ADDITIVE COOLING: start_temp = 0.002, final_temp = 0.00025
  * LINEAR ADDITIVE COOLING: start_temp = 0.003, final_temp = 0.0005
  * LINEAR MULTIPLICATIVE COOLING: start_temp = 0.0025, final_temp = 0.0005

For example, the command used to analyst file **example.phy** by MPBoot-SA with EXPONENTIAL ADDITIVE COOLING scheduler, starting temperature 9.2 and minimum temperature 5.3 is:

`mpboot-avx -s example.phy -pure_sa -cooling_schedule 1 -start_temp 9.2 -final_temp 5.3`


