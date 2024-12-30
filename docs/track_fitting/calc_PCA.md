# PCA Functions in PR3DCluster

The PR3DCluster class implements several Principal Component Analysis (PCA) functions for analyzing geometric data in 3D space. Here's a detailed breakdown of each function:

## 1. `void Calc_PCA()`

This is the main PCA calculation function that operates on the entire cluster.

### Key Features:
- Calculates the center of mass of all points in the cluster
- Constructs a 3x3 covariance matrix from the point cloud
- Uses TMatrixDEigen to compute eigenvalues and eigenvectors
- Stores results in:
  - `center`: Center of mass
  - `PCA_axis[3]`: Three principal axes
  - `PCA_values[3]`: Corresponding eigenvalues

### Implementation Details:
```cpp
// Calculate center of mass
for (auto it = mcells.begin(); it!=mcells.end();it++){
    PointVector ps = (*it)->get_sampling_points();
    for (int k=0;k!=ps.size();k++){
        center.x += ps.at(k).x;
        center.y += ps.at(k).y;
        center.z += ps.at(k).z;
        nsum ++;
    }
}
```

## 2. `void Calc_PCA(PointVector& points)`

An overloaded version that performs PCA on a specific set of points.

### Key Features:
- Takes an input vector of points instead of using the cluster's points
- Similar algorithm to the main Calc_PCA()
- Useful for analyzing subsets of points or temporary point collections

### Key Differences from Main Version:
- Works with provided points rather than cluster's mcells
- More flexible for targeted analysis

## 3. `TVector3 calc_PCA_dir(Point& p, double dis)`

Calculates PCA direction for points within a sphere around a given point.

### Key Features:
- Takes a center point and radius (dis)
- Only considers points within the specified radius
- Returns a TVector3 representing the primary PCA direction
- Uses charge weighting in calculations

### Usage:
```cpp
// Calculate main direction within 5cm radius
TVector3 direction = calc_PCA_dir(some_point, 5.0*units::cm);
```

## 4. `TVector3 calc_PCA_dir(Point& p, PointVector& pts)`

Another overload that calculates PCA direction for a specific set of points relative to a center point.

### Key Features:
- Takes a center point and vector of points
- No distance filtering (uses all provided points)
- Returns primary PCA direction as TVector3
- Useful for pre-filtered or selected point sets

## Common Elements Across Functions

### Covariance Matrix Calculation:
All functions use a similar approach to build the covariance matrix:
```cpp
cov_matrix(0,0) = sum((x - center.x)^2)
cov_matrix(0,1) = sum((x - center.x)(y - center.y))
cov_matrix(0,2) = sum((x - center.x)(z - center.z))
// ... etc for other elements
```

### Eigenvector Processing:
```cpp
// Normalize eigenvectors
for (int i=0;i!=3;i++){
    PCA_axis[i].x = eigen_vectors(0,i)/magnitude;
    PCA_axis[i].y = eigen_vectors(1,i)/magnitude;
    PCA_axis[i].z = eigen_vectors(2,i)/magnitude;
}
```

## Usage Considerations

1. **When to Use Each Version:**
   - `Calc_PCA()`: For analyzing entire clusters
   - `Calc_PCA(PointVector&)`: For custom point sets
   - `calc_PCA_dir(Point&, double)`: For local direction analysis
   - `calc_PCA_dir(Point&, PointVector&)`: For direction with pre-selected points

2. **Performance Considerations:**
   - Main `Calc_PCA()` processes all cluster points - can be expensive for large clusters
   - Versions with point selection are more efficient for local analysis
   - Consider point count when choosing between versions

3. **Error Handling:**
   - All versions check for minimum point count (n≥3)
   - Return zero vectors/values if insufficient points
   - Handle degenerate cases (coplanar/collinear points)