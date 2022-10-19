using DynamicalSystems
using GLMakie
using GeometryBasics
using DrWatson
using CSV
using DataFrames

function emToCentroidsInteractive(fileName)
    mkpath(datadir("exp_pro","emCentroidsInteractive"))

    imageIn = load(fileName)

    fig = GLMakie.Figure()
    ax = GLMakie.Axis(fig[1,1],aspect=DataAspect())
    Makie.deactivate_interaction!(ax, :rectanglezoom)
    hidedecorations!(ax)
    hidespines!(ax)
    image!(ax, rotr90(imageIn))
    GLMakie.display(fig)
    Makie.deactivate_interaction!(ax, :rectanglezoom)

    centroidLocationsObs = Point2[]
    centroidLocationsObs = Observable(centroidLocationsObs) 
    scatter!(ax, centroidLocationsObs, color = :orange, markersize=20)

    spoint = select_point(ax.scene)

    on(spoint) do z
        x, y = z
        push!(centroidLocationsObs[],z)
        centroidLocationsObs[] = centroidLocationsObs[]
    end

    removeButton = Button(fig[2,1][1,1]; label = "remove last", tellwidth = false)
    on(removeButton.clicks) do clicks; pop!(centroidLocationsObs[]); centroidLocationsObs[]=centroidLocationsObs[]; end

    saveButton = Button(fig[2,1][1,2]; label = "save", tellwidth = false)
    on(saveButton.clicks) do clicks; centroidLocations=centroidLocationsObs[]; save(datadir("exp_pro","emCentroidsInteractive","$(splitpath(fileName)[end][1:end-4]).jld2"),@strdict centroidLocations); end
end

runs = Vector(readdlm(datadir("exp_pro","filesToUse.txt"))[:,1])
lengthMeasurements = DataFrame(CSV.File(datadir("exp_pro","lengthMeasurements","lengthMeasurements.csv")))

# for r in runs[1:2]
#     emToCentroidsInteractive(datadir("exp_pro","cropped",r))
# end