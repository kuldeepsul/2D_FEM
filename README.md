##### This solver performs linear static analysis of frame structures based on Euler-Bernoulli Beam Theory, under the assumption of small deformations.

### 1. Governing Equation

$$
\begin{aligned}
Ku = f \\
\end{aligned}
$$

$K = Global Stiffness Matrix$
$u = Displacement Vector$
$f = Force Vector$


### 2. Element Formulation
$$
\mathbf{K}_{E,local} = \frac{E}{L} \begin{bmatrix}
A & 0 & 0 & -A & 0 & 0 \\
0 & \frac{12I}{L^2} & \frac{6I}{L} & 0 & -\frac{12I}{L^2} & \frac{6I}{L} \\
0 & \frac{6I}{L} & 4I & 0 & -\frac{6I}{L} & 2I \\
-A & 0 & 0 & A & 0 & 0 \\
0 & -\frac{12I}{L^2} & -\frac{6I}{L} & 0 & \frac{12I}{L^2} & -\frac{6I}{L} \\
0 & \frac{6I}{L} & 2I & 0 & -\frac{6I}{L} & 4I
\end{bmatrix}
$$

### 3. Co-ordinate Transformation
##### Element Local Co-ordinates
- The stiffness matrix achieved from the element formulation is based on the co-ordinate system that aligns to the lateral and longitudinal directions of element.
##### Model Global Co-ordinates
- The global is based on the global co-ordinate system of the whole model, therefore the stiffness of element assembled in the global stiffness matrix must be in global directions.

![[Pasted image 20260210150429.png]]

##### Transformation Matrix

- A  $6*6$  transformation matrix is used because , each node consists 3 degrees of freedom , and each element contain 2 nodes.

$$
\alpha = Angle\ of\ element to\ the\ global\ x-direction\  ; 
c = \cos \alpha\  ;
s = \sin \alpha 
$$

$$
T = \begin{bmatrix}
c  & -s & 0 & 0 & 0 & 0  \\
s  & c & 0 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 & 0 & 0 &  \\
0 & 0 & 0 & c & -s & 0 \\
0 & 0 & 0 & s & c & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
\end{bmatrix}
$$

##### K-Local to K-Global Transformation

$$
K_{Eglobal} = T^TK_{Elocal}T
$$

$T^T$ = *Transpose of Transformation Matrix*

### 4. Global Stiffness Matrix Assembly

- $K_{Eglobal}$ for each element is assembled into the $K_{Mglobal}$ the global stiffness matrix of the model.
- Each element contains a DOF map containing the DOF id's of the $K_{Eglobal}$ rows and columns.

``` c++
// assume Global dof's for a model are {0,1,2,3,4,5,6,7,8,9,10}.
// K-global is 10*10 matrix.

// and some element connects node02 and node03,where dof for node02 are {3,4,5}
// and {6,7,8}. Since the element has 6 dofs (three at each node) the element    // stiffness matrix is a 6*6 matrix.

// the dof map inside the element object will be a array of global dof ids for the // element. for our case dof_map = {3,4,5,6,7,8}.

			 dof3    dof4    dof5    dof6   dof7   dof8
		/*	 ke11    ke12    ke13    ke14   ke15   ke16    */ dof3
		/*	 ke21    ke22    ke23    ke24   ke25   ke26    */ dof4
		/*	 ke31    ke32    ke33    ke34   ke35   ke36    */ dof5
		/*	 ke41    ke42    ke43    ke44   ke45   ke46    */ dof6
		/*	 ke51    ke52    ke53    ke54   ke55   ke56    */ dof7
		/*	 ke61    ke62    ke63    ke64   ke65   ke66    */ dof8
		
// Global_K will contain a larger stiffness matrix having dof_ids as its indices.
// the above element will populate the 6*6 patch of Global_K i.e. from 3rd row to // 8th row , and from 3rd column to 8th column.

// If some other element contain the dofs same as our previous elements the values 
// populated by this new elements will be added to the values populated by 
// previous elements.

for(element : all elements)
{
	for( i < size(dof_map) )
	{
		for(j < size(dof_map))
		{
			Global_K[dof_map[i],dof_map[j]] += Element_K[i,j];
		}
	}
}
```


### 5. Boundary Conditions

- Boundary are applied by the method of elimination.
- Rows and Columns of the Global stiffness matrix are eliminated based on boundary conditions prescribed in the input file , and a Reduced Global stiffness matrix is formed.
- Ex.
```
*BC
1,1,0          // Dof - 1 of node id - 1 , is fixed .
1,2,0          // Dof - 2 of node id - 1 , is fixed .
1,3,0          // Dof - 3 of node id - 1 , is fixed .
```


### 6. Linear System Solution
- PLU Decomposition is used as linear system solver.
- A dense PLU Decomposition is currently used for clarity and educational purposes; sparse solvers are planned for scalability.

$$
\begin{aligned}
Ku=f \\
K=PLU \\
Pa=f \\
Lb=a \\
Uu=b
\end{aligned}
$$

### 7. Input / Output File Format

#### Input
- Input file format - .inp 
- Input File must be in this specific order . Material -> Nodes -> Frames -> Boundary Conditions  ->  Forces.
- Material -> Material ID, Young's Modulus , Poisson's Ratio
- Node -> Node ID, X , Y
- Frame -> Element ID, Node01 ID , Node02 ID , Area , Moment of Inertia , Material ID
- BC -> Node ID,DOF ID, Value
- Force -> Node ID,DOF ID, Value

Ex. Input for a two element cantilever beam.
```
*Material
1,210000000,0.28
*Node
1,0,0
2,5,0
3,10,0
*Frame
1,1,2,0.1,0.001,1
2,2,3,0.1,0.001,1
*BC
1,1,0
1,2,0
1,3,0
*Force
3,2,-10000
```

#### Output
- Results are output as a .out file.
- Format - Node ID , DOF ID , Value

### 8. Result Verification

#### Test case 1: Pure Axial Tension.
##### Setup -
- Node 1 : ( 0 $m$ , 0 $m$) - Fixed .
- Node 2 : ( 10 $m$ , 0 $m$) - Free.
- Material Properties - E = 200 * $10^9$  $\frac{N} {m^2}$ , $\nu$ = 0.28 .
- Frame Properties - A = 0.01 $m^2$  , I = 0.0001 $m^{4}$.
- Load - 100,000 $N$.
##### Analytical result - 
- Formula ->  $dl$ = $\frac{{F.L}}{A.E}$ = $0.0005m$
- zero lateral displacement and angle.
##### Test result - 
```
1,1,0
1,2,0
1,3,0
2,1,0.005
2,2,0
2,3,0
```

##### Result Plots
![[Pasted image 20260211131022.png]]

#### Test case 2: Pure Bending.

##### Setup -
- Node 1 : ( 0 $m$ , 0 $m$) - Fixed .
- Node 2 : ( 2 $m$ , 0 $m$) - Free.
- Material Properties - E = 200 * $10^9$  $\frac{N} {m^2}$ , $\nu$ = 0.28 .
- Frame Properties - A = 0.01 $m^2$  , I = 0.0001 $m^{4}$.
- Load - -10,000 $N$. (Downwards)
##### Analytical result - 
- Formula 
$v_{max}$ = ${P.L^3} /{3.E.I}$ = -0.0013$m$
$\theta_{max}$ =  ${P.L^2} /{2.E.I}$ = -0.001$rad$
- Zero longitudinal displacement.

##### Test result - 
```
1,1,0
1,2,0
1,3,0
2,1,0
2,2,-0.00133333
2,3,-0.001
```

##### Result plots
![[Pasted image 20260211131240.png]]
Scaled - 100 times
![[Pasted image 20260211131353.png]]

#### Test case 3: Tension on Frame at an Angle.

##### Setup -
- Node 1 : ( 0 $m$ , 0 $m$) - Fixed .
- Node 2 : ( 1 $m$ , 1 $m$) - Free.
- Material Properties - E = 200 * $10^9$  $\frac{N} {m^2}$ , $\nu$ = 0.28 .
- Frame Properties - A = 0.01 $m^2$  , I = 0.0001 $m^{4}$.
- Load :- 10,000 $N$. (Upwards) and 10,000 $N$ . (Forward)

Analytical result - 
- Formula ->  $dl$ = $\frac{{F.L}}{A.E}$ =  0.00001$m$ =>$dl\cos{\alpha}$ = $dx$ = $dy$ =  $7.071*e^{-6}m$
- Longitudinal tension , displacement in both x and y direction must be identical.
- Zero angle of rotation. 

##### Test result - 
```
1,1,0
1,2,0
1,3,0
2,1,7.07107e-06
2,2,7.07107e-06
2,3,0
```

##### Result Plot
![[Pasted image 20260211144409.png]]

### 9. Software Architecture

- `Model.cpp`  - Classes  Node, Element ,Model.
- `Material.cpp` - Contains a very minimal Material implementation.
- `Matrix.cpp` -  Matrix operations and linear system solver.
- `Input_Reader.cpp` - functions to input and output data.
- `Main.cpp` - Main loop.

### 10. Limitation

- No geometric non- linearity.
- No material non-linearity.
- No dynamics.

### 11. Roadmap

- Iso- Parametric Quad implementation.
- Material Non-Linearity.

### Overview
->Read Input.
->Generate Model.
->Generate Element Stiffness Matrix.
->Transform Element Stiffness Matrix to Elements oriented at arbitrary angle.
->Assemble Global Stiffness Matrix.
->Apply Boundary Condition.
->Solve linear system [K]{x} = {f} , by PLU Decomposition.
->Output Results.

