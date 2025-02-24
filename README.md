## 2D Truss Analyzer
2D Truss Anlalyzer is a student-level application written in C++ featuring DXF output generator for initial and deformed shape of the truss.

### Features

- Based on Matrix Analysis
- Unloaded vs loaded truss viewer
- Standard DXF format output

![Truss Analysis Result: Deflection](https://github.com/user-attachments/assets/362af45c-06db-4bf9-ba86-3e1e39399371)

### Input File Format
The application reads inputes form a `TA.input` text file in separate lines. The input format described below.

####File Format
```javascript
8                                   // Truss Color Code
One_Color                           // Truss Display Option
3                                   // Number of Joints
0                                   // Joint1 X
0                                   // Joint1 Y
0                                   // Joint2 X
300                                 // Joint2 Y
400                                 // Joint3 X
300                                 // Joint3 Y
3                                   // Number of Members
2400                                // Members' limit stress: fy
1                                   // Joint A ID
2                                   // Joint B ID
3                                   // Member's Area
2000000                             // Member's Module of Elasticity
2                                   // Joint A ID
3                                   // Joint B ID
4                                   // Member's Area
2000000                             // Member's Module of Elasticity
3                                   // Joint A ID
1                                   // Joint B ID
5                                   // Member's Area
2000000                             // Member's Module of Elasticity
2                                   // Number of Supports
1                                   // Joint ID
1                                   // Horizontal Displacement
1                                   // Vertical Displacement
2                                   // Joint ID
1                                   // Horizontal Displacement
1                                   // Vertical Displacement
1                                   // Number of Loads
3                                   // Joint ID
707.1                               // Horizontal Force
-707.1                              // Vertical Force
```
