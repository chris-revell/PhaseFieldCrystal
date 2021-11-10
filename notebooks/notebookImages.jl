### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 30895b10-cd1b-11eb-13f1-e75fc2519aa0
using Images, ImageBinarization,PlutoUI

# ╔═╡ d375a2f1-133a-4d08-96ed-b4bd1d446dce
function myBinarise(image)
	if image.val<0.8
		return Gray(0.0)
	else
		return Gray(1.0)
	end
end

# ╔═╡ e3f67fa5-c972-4381-993d-3929afeb50cb
@bind distance Slider(1:50)

# ╔═╡ 244b65d1-2e50-460c-b3e4-04cc927e5da7
begin
	#image = load("/Users/christopher/Dropbox (The University of Manchester)/Other/EM-images/e13.5tail(TS)/635.tif")[:,1:4500]
	image = Gray.(load("/Users/christopher/Desktop/635.jpg"))
	#eroded = erode(image)
	#dilated = dilate(eroded)
	gaussianFiltered = imfilter(image,Kernel.gaussian(distance))
	myBinarised = myBinarise.(gaussianFiltered)#,MinimumError())
end

# ╔═╡ bf127322-fb08-4665-89bc-3e7d798bae48


# ╔═╡ 15f81d45-3b20-4526-927f-23d103ba7ce7


# ╔═╡ ace0008d-d36a-4601-8c2f-574b7b830fa0
#myBinarised = binarize(image,Otsu())

# ╔═╡ Cell order:
# ╠═30895b10-cd1b-11eb-13f1-e75fc2519aa0
# ╠═d375a2f1-133a-4d08-96ed-b4bd1d446dce
# ╠═e3f67fa5-c972-4381-993d-3929afeb50cb
# ╠═244b65d1-2e50-460c-b3e4-04cc927e5da7
# ╠═bf127322-fb08-4665-89bc-3e7d798bae48
# ╠═15f81d45-3b20-4526-927f-23d103ba7ce7
# ╠═ace0008d-d36a-4601-8c2f-574b7b830fa0
