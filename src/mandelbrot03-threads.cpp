#include <iostream>
#include <string>
#include <cxxopts.hpp>
#include <thread>
#include <atomic>

std::string title("Mandelbrot_threads");

#include "common.inc"

double mylog2(double value)
{
    constexpr int mantissa_bits = 52, exponent_bias = 1022;
    const double  half         = 0.5;
    std::uint64_t half_bits    = reinterpret_cast<const std::uint64_t&>(half);
    int e,lt;
    uint64_t m;
    double x, dbl_e, z, y, u, t;
    m = reinterpret_cast<const std::uint64_t&>(value);
    e = m >> mantissa_bits; // frexp(). e = exponent, m = mantissa
    m &= std::uint64_t((1ull << mantissa_bits)-1);
    m |= half_bits;
    x = reinterpret_cast<const double&>(m);
    lt = (x < 1/std::sqrt(2.)) ? -1 : 0;
    dbl_e = e + lt - exponent_bias;
    z = x - (half + (lt ? 0. : half));
    y = half * (x - (lt ? half : 0.)) + half;
    x = z/y;
    z = x*x;
    u = z   + -3.56722798512324312549E1;
    t =       -7.89580278884799154124E-1;
    u = u*z +  3.12093766372244180303E2;
    t = t*z +  1.63866645699558079767E1;
    u = u*z + -7.69691943550460008604E2;
    t = t*z + -6.41409952958715622951E1;
    y = z* (t/u) + (half+half);
    return x*(y*std::log2(std::exp(1.))) + dbl_e;
}

template<bool WithMoment>
double Iterate(double zr, double zi)
{
    const double escape_radius_squared = ESCAPE_RADIUS_SQUARED;
    const int maxiter = MAXITER;
    double cr = zr, sr = cr;
    double ci = zi, si = ci;
    double dist;
    int iter = maxiter, notescaped = -1;

    if(zr*(1+zr*(8*zr*zr+(16*zi*zi-3)))+zi*zi*(8*zi*zi-3) < 3./32 || ((zr+1)*(zr+1)+zi*zi)<1./16) { iter=0; }

    while(notescaped)
    {
        double r2 = cr * cr;
        double i2 = ci * ci;
        dist = r2 + i2;

        notescaped &= ((iter != 0) & (dist < escape_radius_squared)) ? -1 : 0;
        iter += notescaped;

        double ri = cr * ci;
        ci = zi + (ri * 2);
        cr = zr + (r2 - i2);

        if(WithMoment)
        {
            bool notmoment = iter & (iter-1);
            iter = (cr == sr && ci == si) ? 0 : iter;
            sr = notmoment ? sr : cr;
            si = notmoment ? si : ci;
        }
    }
    return iter ? mylog2( maxiter-iter + 1 - mylog2(mylog2(dist) / 2)) * (4/std::log2(std::exp(1.))) : 0;
}
int main(int argc, char** argv) {
    bool timings = false;
    int retval = 0;
    if(!parseOptions(argc, argv, title, retval, timings))
        return retval;
	
    std::vector<float> pixels(width * height * 3);

    Window window(width, height, title, timings);
	bool needMoment = true;

    std::cout << "Hardware threads used: " << std::thread::hardware_concurrency() << "\n";

    while(!window.shouldClose()) {

        double zr, zi, xscale, yscale; 

        // Set coordinates and scale:
        if(timings)
            setCoordinatesFrame(zr, zi, xscale, yscale, width, height, window.getFrameCount());
        else
            setCoordinatesTime(zr, zi, xscale, yscale, width, height);
	
        std::atomic<unsigned>    y_done{0}, n_inside{0};
        std::vector<std::thread> threads;
        
        for(unsigned n=0; n<std::thread::hardware_concurrency(); ++n)
			threads.emplace_back([&](){
                unsigned count_inside = 0;
                for(unsigned y; (y = y_done++) < height; )
                {
                    double i = zi+yscale*int(y-height/2);
                    if(needMoment)
                        for(unsigned x=0; x<width; ++x)
                        {
                            double v = Iterate<true>( zr+xscale*int(x-width/2), i );
                            if(v == 0.) ++count_inside;
                            auto index = x*3 + y*width*3;
							auto color = getColor(x, y, v, false);
                            pixels[index] = std::get<0>(color); //x;
                    		pixels[index + 1 ] = std::get<1>(color); //= y;
                    		pixels[index + 2] = std::get<2>(color);
                        }
                    else
                        for(unsigned x=0; x<width; ++x)
                        {
                            double v = Iterate<false>( zr+xscale*int(x-width/2), i );
                            if(v == 0.) ++count_inside;
                        	auto index = x*3 + y*width*3;
                            auto color = getColor(x, y, v, false);
                            pixels[index] = std::get<0>(color); //x;
                            pixels[index + 1 ] = std::get<1>(color); //= y;
                            pixels[index + 2] = std::get<2>(color);
						}
                }
                n_inside += count_inside;
            });

        for(auto& t: threads) t.join();
    
        needMoment = n_inside >= (width*height)/1024;
        window.draw(pixels, width, height);

    }
    

    return 0;
}
