function defineInternalCoordinates(zMatrixCoordinates::Vector{Float64})::Vector{Float64}
    convertToRadians::Float64 = pi/180
    # Trivial case - no dihedrals to convert just degrees to radians
    internalCoordinates::Vector{Float64} = zeros(length(zMatrixCoordinates))
    internalCoordinates[1:2] = zMatrixCoordinates[1:2]
    internalCoordinates[3] = zMatrixCoordinates[3]*convertToRadians
    return internalCoordinates
end

function defineXiCoordinates(internalCoordinates::Vector{Float64}, structuralParameters::Vector{Float64})::Vector{Float64}
    convertToRadians::Float64 = pi/180
    numberOfModes::Int64 = length(internalCoordinates)
    numberOfAtoms::Int64 = (numberOfModes + 6)/3
    numberOfStretches::Int64 = numberOfAtoms - 1
    xi::Vector{Float64} = zeros(numberOfModes)
    equilibriumParameters::Vector{Float64} = structuralParameters[1:numberOfModes]
    morseParameters::Vector{Float64} = structuralParameters[numberOfModes+1:end]
    
    # Define Morse-like coordinates for stretches
    stretchDisplacement::Vector{Float64} = internalCoordinates[1:numberOfStretches] - equilibriumParameters[1:numberOfStretches]
    xi[1:numberOfStretches] = 1 .- exp.(-morseParameters.*stretchDisplacement)
    # xi[1:numberOfStretches] = stretchDisplacement

    # xi[numberOfStretches+1:end] = (internalCoordinates[numberOfStretches+1:end] - equilibriumParameters[numberOfStretches+1:end]).*convertToRadians
    # # Trigonometric-type function for the bending 
    xi[numberOfStretches+1:end] = cos.(internalCoordinates[numberOfStretches+1:end].*convertToRadians) .- cos.(equilibriumParameters[numberOfStretches+1:end].*convertToRadians)
    # xi[numberOfStretches+1:end] = cos.(pi .- internalCoordinates[numberOfStretches+1:end].*convertToRadians) .- cos.(pi .- equilibriumParameters[numberOfStretches+1:end].*convertToRadians)
    # xi[numberOfStretches+1:end] = sin.(pi .- internalCoordinates[numberOfStretches+1:end].*convertToRadians) .- cos.(pi .- equilibriumParameters[numberOfStretches+1:end].*convertToRadians)
    return xi
end