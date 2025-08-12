public class FluidSimulation
{

    public enum FieldType
    {
        X_VELOCITY_FIELD = 0,
        Y_VELOCITY_FIELD = 1,
        SMOKE_FIELD = 2
    }
    // Grid Dimensions
    private  int gridWidth;
    private  int gridHeight;
    private  int numCells;
    private  float cellSize;
    private  float inverseCellSize;
    private  float halfCellSize;

    // Velocity Fields
    private float[] velocityX;
    private float[] velocityY;
    private float[] newVelocityX;
    private float[] newVelocityY;

    // Scalar Fields
    private float[] pressure;
    private float[] smokeDensity;
    private float[] newSmokeDensity;
    private float[] cellType;

    // Simulation Parameters
    private float density;

    // Helper method for array indexing
    private int GridIndex(int gridX, int gridY) => gridX * gridHeight + gridY;

    // Methods for Rendering
    public float[] Pressure => pressure;
    public float[] SmokeDensity => smokeDensity;
    public float[] VelocityX => velocityX;
    public float[] VelocityY => velocityY;
    public float[] CellType => cellType;

    // Grid properties for rendering
    public int GridWidth => gridWidth;
    public int GridHeight => gridHeight;
    public float CellSize => cellSize;

    public void setCellSize(float newCellSize)
    {
        this.cellSize = newCellSize;
        this.inverseCellSize = 1.0f / newCellSize;
        this.halfCellSize = 0.5f * newCellSize;
    }

    public void setGridSize(int newGridWidth, int newGridHeight)
    {
        this.gridWidth = newGridWidth + 2; // Adding ghost cells
        this.gridHeight = newGridHeight + 2; // Adding ghost cells
        this.numCells = this.gridWidth * this.gridHeight;

        // Reinitialize arrays with the new grid size
        this.velocityX = new float[this.numCells];
        this.velocityY = new float[this.numCells];
        this.newVelocityX = new float[this.numCells];
        this.newVelocityY = new float[this.numCells];
        this.pressure = new float[this.numCells];
        this.smokeDensity = new float[this.numCells];
        this.newSmokeDensity = new float[this.numCells];
        this.cellType = new float[this.numCells];

        for (int i = 0; i < this.numCells; i++)
        {
            this.smokeDensity[i] = 1.0f; // Initialize all cells as fluid
        }
    }

    // CellType determines if a ce|
    public void setCellType(int gridX, int gridY, float type)
    {
        if (gridX < 0 || gridX >= gridWidth || gridY < 0 || gridY >= gridHeight)
            throw new ArgumentOutOfRangeException("Grid coordinates are out of bounds.");

        int index = GridIndex(gridX, gridY);
        this.cellType[index] = type;
    }

    public FluidSimulation(float density, int gridWidth, int gridHeight, float cellSize)
    {
        this.density = density;
        this.gridWidth = gridWidth + 2; // Adding ghost cells
        this.gridHeight = gridHeight + 2; // Adding ghost cells
        this.numCells = this.gridWidth * this.gridHeight;
        this.cellSize = cellSize;
        this.inverseCellSize = 1.0f / cellSize;
        this.halfCellSize = 0.5f * cellSize;

        this.velocityX = new float[this.numCells];
        this.velocityY = new float[this.numCells];
        this.newVelocityX = new float[this.numCells];
        this.newVelocityY = new float[this.numCells];
        this.pressure = new float[this.numCells];
        this.smokeDensity = new float[this.numCells];
        this.newSmokeDensity = new float[this.numCells];
        this.cellType = new float[this.numCells];
        for (int i = 0; i < this.numCells; i++)
        {
            this.smokeDensity[i] = 1.0f; // Initialize all cells as fluid
        }
    }

    public void Integrate(float deltaTime, float gravity)
    {
        for (int gridX = 1; gridX < gridWidth; gridX++)
        {
            for (int gridY = 1; gridY < gridHeight - 1; gridY++)
            {
                int currentCell = GridIndex(gridX, gridY);
                int cellBelow = GridIndex(gridX, gridY - 1);


                if ((cellType[currentCell] != 0.0f) && (cellType[cellBelow] != 0.0f))
                {
                    velocityY[currentCell] += gravity * deltaTime;
                }
            }
        }
    }

    public void SolveIncompressibility(int numIterations, float deltaTime, float overRelaxation)
    {
        float pressureCoefficient = this.density * this.cellSize / deltaTime;

        for (int iter = 0; iter < numIterations; iter++)
        {
            for (int i = 1; i < gridWidth - 1; i++)
            {
                for (int j = 1; j < gridHeight - 1; j++)
                {
                    int currentCell = GridIndex(i, j);
                    if (cellType[currentCell] == 0.0f)
                        continue;

                    float neighborSum = cellType[currentCell];
                    float leftNeighbor = cellType[GridIndex(i - 1, j)];
                    float rightNeighbor = cellType[GridIndex(i + 1, j)];
                    float belowNeighbor = cellType[GridIndex(i, j - 1)];
                    float aboveNeighbor = cellType[GridIndex(i, j + 1)];
                    neighborSum = leftNeighbor + rightNeighbor + belowNeighbor + aboveNeighbor;
                    if (neighborSum == 0.0f)
                        continue;

                    float divergence = velocityX[GridIndex(i + 1, j)] - velocityX[currentCell] +
                                velocityY[GridIndex(i, j + 1)] - velocityY[currentCell];

                    float cellPressure = -divergence / neighborSum;
                    cellPressure *= overRelaxation;
                    this.pressure[currentCell] += pressureCoefficient * cellPressure;

                    velocityX[currentCell] -= leftNeighbor * cellPressure;
                    velocityX[GridIndex(i + 1, j)] += rightNeighbor * cellPressure;
                    velocityY[currentCell] -= belowNeighbor * cellPressure;
                    velocityY[GridIndex(i, j + 1)] += aboveNeighbor * cellPressure;
                }
            }
        }
    }

    // Extrapolate boundry velocities to ensure they are consistent with the surrounding cells
    // This is useful for maintaining fluid flow at the edges of the grid
    public void Extrapolate()
    {
        for (int i = 0; i < gridWidth; i++)
        {
            velocityX[GridIndex(i, 0)] = velocityX[GridIndex(i, 1)];
            velocityX[GridIndex(i, gridHeight - 1)] = velocityX[GridIndex(i, gridHeight - 2)];
        }
        for (int j = 0; j < gridHeight; j++)
        {
            velocityY[GridIndex(0, j)] = velocityY[GridIndex(1, j)];
            velocityY[GridIndex(gridWidth - 1, j)] = velocityY[GridIndex(gridWidth - 2, j)];
        }
    }

    // Bilinear interpolation to sample the field at a given point
    public float SampleField(float gridX, float gridY, FieldType fieldType)
    {

        // Clamp gridX and gridY to the valid range
        gridX = Math.Max(Math.Min(gridX, this.gridWidth * this.cellSize), this.cellSize);
        gridY = Math.Max(Math.Min(gridY, this.gridHeight * this.cellSize), this.cellSize);

        // Initialize dx and dy (offsets from grid coordinates) based on the field type
        float dx = 0.0f;
        float dy = 0.0f;
        float[] field;

        switch (fieldType)
        {
            case FieldType.X_VELOCITY_FIELD: // U_FIELD
                field = this.velocityX;
                dy = halfCellSize;
                break;
            case FieldType.Y_VELOCITY_FIELD: // V_FIELD
                field = this.velocityY;
                dx = halfCellSize;
                break;
            case FieldType.SMOKE_FIELD: // S_FIELD
                field = this.smokeDensity;
                dx = halfCellSize;
                dy = halfCellSize;
                break;
            default:
                throw new ArgumentException("Invalid field type");
        }

        // Find the bounding grid cell indices
        float leftCell = Math.Min((float)Math.Floor((gridX - dx) * inverseCellSize), this.gridWidth - 1);
        float interpolationX = ((gridX - dx) - leftCell * this.cellSize) * inverseCellSize;
        float rightCell = Math.Min(leftCell + 1, this.gridWidth - 1);

        float bottomCell = Math.Min((float)Math.Floor((gridY - dy) * inverseCellSize), this.gridHeight - 1);
        float interpolationY = ((gridY - dy) - bottomCell * this.cellSize) * inverseCellSize;
        float topCell = Math.Min(bottomCell + 1, this.gridHeight - 1);

        // complementary interpolation factors
        float compInterpolationX = 1.0f - interpolationX;
        float compInterpolationY = 1.0f - interpolationY;

        // Bilinear interpolation using the four surrounding grid points
        return
            compInterpolationX * compInterpolationY * field[GridIndex((int)leftCell, (int)bottomCell)] +
            interpolationX * compInterpolationY * field[GridIndex((int)rightCell, (int)bottomCell)] +
            interpolationX * interpolationY * field[GridIndex((int)rightCell, (int)topCell)] +
            compInterpolationX * interpolationY * field[GridIndex((int)leftCell, (int)topCell)];
    }

    // Average velocity in the x-direction at the given grid cell
    public float AverageVelocityX(int gridX, int gridY)
    {
        return (this.velocityX[GridIndex(gridX, gridY - 1)] + this.velocityX[GridIndex(gridX, gridY)] +
                this.velocityX[GridIndex(gridX + 1, gridY - 1)] + this.velocityX[GridIndex(gridX + 1, gridY)]) * 0.25f;
    }

    // Average velocity in the y-direction at the given grid cell
    public float AverageVelocityY(int gridX, int gridY)
    {
        return (this.velocityY[GridIndex(gridX - 1, gridY)] + this.velocityY[GridIndex(gridX, gridY)] +
                this.velocityY[GridIndex(gridX - 1, gridY + 1)] + this.velocityY[GridIndex(gridX, gridY + 1)]) * 0.25f;
    }

    // Advect the velocity fields based on the current velocities
    public void AdvectVelocity(float deltaTime)
    {
        // Copy current velocities to new arrays for updating
        Array.Copy(this.velocityX, this.newVelocityX, this.numCells);
        Array.Copy(this.velocityY, this.newVelocityY, this.numCells);

        // Loop through the grid, excluding the ghost cells (boundary cells)
        for (int gridX = 1; gridX < this.gridWidth - 1; gridX++)
        {
            for (int gridY = 1; gridY < this.gridHeight - 1; gridY++)
            {
                int currentCell = GridIndex(gridX, gridY);

                // x velocity component
                if (this.cellType[currentCell] != 0.0f && this.cellType[GridIndex(gridX - 1, gridY)] != 0.0f && gridY < this.gridHeight - 1)
                {
                    float cellEdge = gridX * cellSize;
                    float cellCenter = gridY * cellSize + halfCellSize;
                    float cellVelocityX = this.velocityX[currentCell];
                    float cellVelocityY = AverageVelocityY(gridX, gridY);
                    cellEdge -= deltaTime * cellVelocityX;
                    cellCenter -= deltaTime * cellVelocityY;
                    cellVelocityX = SampleField(cellEdge, cellCenter, FieldType.X_VELOCITY_FIELD);
                    this.newVelocityX[currentCell] = cellVelocityX;
                }
                // y velocity component
                if (this.cellType[currentCell] != 0.0f && this.cellType[GridIndex(gridX, gridY - 1)] != 0.0f && gridX < this.gridWidth - 1)
                {
                    float cellEdge = gridX * cellSize + halfCellSize;
                    float cellCenter = gridY * cellSize;
                    float cellVelocityX = AverageVelocityX(gridX, gridY);
                    float cellVelocityY = this.velocityY[currentCell];
                    cellEdge -= deltaTime * cellVelocityX;
                    cellCenter -= deltaTime * cellVelocityY;
                    cellVelocityY = SampleField(cellEdge, cellCenter, FieldType.Y_VELOCITY_FIELD);
                    this.newVelocityY[currentCell] = cellVelocityY;
                }
            }
        }
        // Update the velocity fields with the newly computed values
        Array.Copy(this.newVelocityX, this.velocityX, this.numCells);
        Array.Copy(this.newVelocityY, this.velocityY, this.numCells);
    }

    // Advect the smoke density field based on the current velocities
    public void AdvectSmoke(float deltaTime)
    {
        // Copy current smoke density to new array for updating
        Array.Copy(this.smokeDensity, this.newSmokeDensity, this.numCells);


        for (int gridX = 1; gridX < this.gridWidth - 1; gridX++)
        {
            for (int gridY = 1; gridY < this.gridHeight - 1; gridY++)
            {
                if (this.cellType[GridIndex(gridX, gridY)] != 0.0f)
                {
                    float averageVelX = (this.velocityX[GridIndex(gridX, gridY)] + this.velocityX[GridIndex(gridX + 1, gridY)]) * 0.5f;
                    float averageVelY = (this.velocityY[GridIndex(gridX, gridY)] + this.velocityY[GridIndex(gridX, gridY + 1)]) * 0.5f;
                    float tracebackX = gridX * cellSize + halfCellSize - deltaTime * averageVelX;
                    float tracebackY = gridY * cellSize + halfCellSize - deltaTime * averageVelY;

                    this.newSmokeDensity[GridIndex(gridX, gridY)] = SampleField(tracebackX, tracebackY, FieldType.SMOKE_FIELD);
                }
            }
        }
        // Update the smoke density field with the newly computed values
        Array.Copy(this.newSmokeDensity, this.smokeDensity, this.numCells);
    }

    // Main simulation method that integrates forces (gravity), solves incompressibility, extrapolates boundary velocities,
    // and advects both velocity and smoke density fields
    public void Simulate(float deltaTime, int numIterations, float gravity, float overRelaxation)
    {
        Integrate(deltaTime, gravity);
        SolveIncompressibility(numIterations, deltaTime, overRelaxation);
        Extrapolate();
        AdvectVelocity(deltaTime);
        AdvectSmoke(deltaTime);
    }
}


