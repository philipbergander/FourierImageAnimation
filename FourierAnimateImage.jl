using FileIO, Images
using Plots

pngFileName = "M:\\Matlab\\imgArrow.png"
savePath = "C:\\Users\\phiber\\Downloads\\imgArrowAnim.gif"

function main(pngFile, savePath)

    #### Reading of image to create the function in the complex plane ####

    # Load PNG file and convert from RBG to 3-dimensional matrix
    img_RGB = load(File{format"PNG"}(pngFile))
    img_Mat = permutedims(channelview(img_RGB),(2,3,1))

    # pxl is a boolean matrix where the elements are true where the drawing exists
    pxl = img_Mat[:,:,1] .!= 1
    findPxl = findall(pxl)

    # f is an nx2-dimensional matrix where each row is a coordinate for the consecutive step of the animation
    # f can thus be seen as the discrete complex function in the time domain
    f = Array{Int64}(undef,length(findPxl),2) 

    #row (r) and column (c) index of an arbitrary point in the drawing, becomes the starting position of the animation
    r,c = Tuple(findPxl[1]) 
    f[1,:] = [r,c]
    pxl[findPxl[1]] = false #Set the "used" pixel to false to prevent it to be used again

    nearIdx = -1:1

    j = 2
    while true
        # Look at the pixels adjacent to the previous iteration's pixel and find one that is "drawn"
        near = pxl[f[j-1,1] .+ nearIdx, f[j-1,2] .+ nearIdx]
        firstNear = findfirst(near)
        if isnothing(firstNear) # if there is no more pixels, break the loop
            break
        end

        #Set the next coordinate
        rt,ct = Tuple(firstNear)
        r = r+nearIdx[rt]
        c = c+nearIdx[ct]
        f[j,:] = [r,c]
        pxl[r,c] = false #Set the "used" pixel to false to prevent it to be used again
        j += 1
    end

    f = f[begin:j-1,:]

    #Set the center of the png as the origin
    mid = ceil.(size(pxl) ./2)

    f[:,1] = f[:,1] .- mid[1]
    f[:,2] = f[:,2] .- mid[2]

    #### Complex Fourier Series  ####
    t = range(0,1,length = size(f,1)+1) # Discrete time steps
    t = t[begin:end-1]

    N = 200 #Number of non-static coefficients. Shall be even
    cp = Vector{Complex{Float64}}(undef,Int64(N/2)) #Positive coefficients
    cn = Vector{Complex{Float64}}(undef,Int64(N/2)) #Negative coefficients
    c0 = sum(f[:,1] .+ 1im.* f[:,2]) ./ size(f,1) #static coefficient
    
    # Calculation of Fourier coefficients
    for k=eachindex(cp)
        cp[k] = sum((f[:,1] .+ 1im.*f[:,2]) .* exp.(-k*2im*pi.*t))./size(f,1)
        cn[k] = sum((f[:,1] .+ 1im.*f[:,2]) .* exp.(k*2im*pi.*t))./size(f,1)
    end
    
    # Calculation of function in frequency domain
    function vecsumFunc(tt)
        #Function to calculate the sum of the phasors
        F = c0;
        for k = eachindex(cp)
            F = F + cp[k]*exp(k*2im*pi*tt) + cn[k]*exp(-k*2im*pi*tt)
        end
        return F
    end

    #### Plotting ####
    # x and y are the real and imaginary parts, respectively 
    x = real(vecsumFunc.(t))
    y = imag(vecsumFunc.(t))

    println("Animating...")

    anim = @animate for i in 1:length(t)
        fourierplot(x,y,i)

        plot!([0,real(c0)],[0,imag(c0)], arrow=false,linewidth=1,color=:black, label=false)
        vecSum = c0;
        for k = eachindex(cp)
            plot!([real(vecSum),real(vecSum + cp[k]*exp(k*2im*pi*t[i]))],
             [imag(vecSum),imag(vecSum + cp[k]*exp(k*2im*pi*t[i]))],
             arrow=false,linewidth=1,color=:dodgerblue4, label=false)
            vecSum = vecSum + cp[k]*exp(k*2im*pi*t[i])
            
            plot!([real(vecSum),real(vecSum + cn[k]*exp(-k*2im*pi*t[i]))],
             [imag(vecSum),imag(vecSum + cn[k]*exp(-k*2im*pi*t[i]))],
             arrow=false,linewidth=1,color=:darkred, label=false)
            vecSum = vecSum + cn[k]*exp(-k*2im*pi*t[i])
        end
    end

    println("Create gif...")

    gif(anim, savePath, fps = 50)

    return
end

@userplot FourierPlot
@recipe function func(cp::FourierPlot)
    x,y,i = cp.args
    n = length(x)
    inds = circshift(1:n,1-i)
    linewidth --> [range(0,0,length = Int64(floor(n/2))); range(0,10,length = Int64(ceil(n/2)))]
    seriesalpha --> [range(0, 0, length = Int64(floor(n/2))); range(0, 1, length = Int64(ceil(n/2)))]
    aspec_ratio --> 1
    label --> false
    size --> (1080,1080)
    framestyle --> :none
    x[inds], y[inds]
end

main(pngFileName,savePath)
