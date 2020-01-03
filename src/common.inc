#include <chrono>

#include <mork/ui/GlfwWindow.h>
#include <mork/render/Texture.h>
#include <mork/render/Framebuffer.h>
#include <mork/render/VertexBuffer.h>
#include <mork/render/VertexArrayObject.h>

constexpr unsigned int MAXITER=8000;
constexpr unsigned int ESCAPE_RADIUS_SQUARED=6*6;

constexpr unsigned int width = 1024;
constexpr unsigned int height = 800;

double frameTime = 0.0;

const char *vertexShaderSource =
    "layout (location = 0) in vec3 aPos;\n"
    "out vec2 texCoord;\n"
    "void main()\n"
    "{\n"
    "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
    "   texCoord = vec2(aPos.x*0.5 + .5, aPos.y*0.5+0.5);\n"
    "}\0";

const char *fragmentShaderSource =
    "out vec4 FragColor;\n"
    "in vec2 texCoord;\n"
    "\n"
    "uniform sampler2D ourTexture;\n"
	"const uint r[12] = uint[12](0x00u,0x40u,0x7Eu,0x11u,0x16u,0x38u,0xFCu,0xD0u,0x5Fu,0xDCu,0xFFu,0x6Bu);\n"
	"const uint g[12] = uint[12](0x00u,0x40u,0x9Fu,0x90u,0x68u,0xCFu,0xFFu,0x99u,0x00u,0x37u,0x8Eu,0x14u);\n"
	"const uint b[12] = uint[12](0x00u,0xE0u,0xFFu,0x9Fu,0x18u,0x3Fu,0x00u,0x24u,0x09u,0x0Au,0xFEu,0xBCu);\n"
	"const int k = 1;\n"
	"const int m = 0x3F;\n"    
"\n"
    "void main()\n"
    "{\n"
	"   FragColor = texture(ourTexture, texCoord);\n"
 	"   int x = int(FragColor.x);\n"
	" 	int y = int(FragColor.y);\n"
	"   float iter = FragColor.z;\n"
	"   //double d = ((((x&4)/4u + (x&2)*2u + (x&1)*16u) + (((x^y)&4)/2u + ((x^y)&2)*4u + ((x^y)&1)*32u))&m)/64.;\n"
    "   float R = mix(r[int(iter)%12], r[int(iter+1)%12], iter-int(iter))/255.0;\n"
	"   float G = mix(g[int(iter)%12], g[int(iter+1)%12], iter-int(iter))/255.0;\n"
	"   float B = mix(b[int(iter)%12], b[int(iter+1)%12], iter-int(iter))/255.0;\n"
	"   FragColor = vec4(R, G, B, 1.0);\n"
	"}\n\0";
                                                                                    
bool parseOptions(int argc, char** argv, const std::string& title, int& retval, bool& timings) {     
    cxxopts::Options options(argv[0], title);
    options.add_options()
        ("h,help", "Print help")
        ("v,verbose", "Verbose processing")
        ("t,timing", "process all frames and compute timings (default false)", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
		("positional",
       "Positional arguments: these are the arguments that are entered "
       "without an option", cxxopts::value<std::vector<std::string>>());
        ;
 	auto result = options.parse(argc, argv);
	  

    timings = result["timing"].as<bool>(); 

    if(result.count("help"))
    {
        std::cout << options.help({""}) << std::endl;
        retval = 0;
        return false;
    }
 
    return true;
}

double GetTime()
{
    static std::chrono::time_point<std::chrono::system_clock> begin = std::chrono::system_clock::now();
    return std::chrono::duration<double>( std::chrono::system_clock::now() - begin ).count();
}


// an opengl window with external render loop
class Window : public mork::GlfwWindow {
    public:
    Window(unsigned int width, unsigned int height, const std::string& title, bool timings) : GlfwWindow(mork::Window::Parameters().size(width, height).name(title)), tex(width, height, GL_RGBA, false), prog(vertexShaderSource, fragmentShaderSource), timings(timings)

    {
        mork::GlfwWindow::waitForVSync(false); 
    
        this->reshape(this->getWidth(), this->getHeight());

        std::vector<mork::vertex_pos3> vertices = {
            mork::vec3f(1.0f,  1.0f, 0.0f),  // top right
            mork::vec3f(1.0f, -1.0f, 0.0f),  // bottom right
            mork::vec3f(-1.0f, -1.0f, 0.0f),  // bottom left
            mork::vec3f(-1.0f,  1.0f, 0.0f)   // top left
        };

        std::vector<unsigned int> indices = {
            0, 1, 3,  // first Triangle
            1, 2, 3   // second Triangle
        };

        vao.bind();
        vb.bind();
        vb.setData(vertices);
        vb.setAttributes();

        ib.bind();
        ib.setData(indices);
        vb.unbind();
        vao.unbind();

        if(timings) {
            std::ofstream of(title + ".timings");
            of << "Frame_no,frametime\n";
            of.close();

        }


        frameTime = GetTime();
    }

    void    draw(const std::vector<float>& pixels, unsigned int width, unsigned int height) {
    
        auto currTime = GetTime();
        frameTime = currTime - frameTime;
        
        mork::Framebuffer::getDefault().bind();
        mork::Framebuffer::getDefault().clear();



        // Copy pixels to GPU:
        tex.bind(0);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGB, GL_FLOAT, (void*)pixels.data());
        
        prog.use();
        prog.getUniform("ourTexture").set(0);
        vao.bind();
        tex.bind(0);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


        GlfwWindow::redisplay(t, dt);

        glfwPollEvents();

        if(timings) {
            std::ofstream of(title + ".timings", std::ios_base::app);
            of << this->getFrameCount() << "," << frameTime << "\n";
            of.close();

        }

        frameTime = GetTime();

    }

    bool shouldClose() const {
        if(this->getFrameCount()>3600)
            return true;
        else
            return glfwWindowShouldClose((GLFWwindow*)glfwWindowHandle)!=0;
    }

    virtual void fpsUpdatedEvent() {
        std::stringstream os;
        os << this->getTitle() + " (FPS: " << this->getFps() << ", frame: " << this->getFrameCount() << ", frameTime: " << frameTime << ")";
        this->setWindowTitle(os.str());
    }

    virtual void reshape(int x, int y)
    {
        // Resizes default framebuffer (will also sort out glViewport
        mork::Framebuffer::getDefault().setSize(mork::vec2i(x, y));
        GlfwWindow::reshape(x, y);
        //tex = mork::Texture<2>(x, y, GL_RGBA, false);
        idle(false);
                
    }
    
    private:
        mork::Texture<2> tex;
        mork::VertexBuffer<mork::vertex_pos3> vb;
        mork::VertexArrayObject vao;
        mork::GPUBuffer<unsigned int, GL_ELEMENT_ARRAY_BUFFER, GL_STATIC_DRAW> ib;

        mork::Program prog;

        bool timings;
};


void setCoordinatesFrame(double& zr, double& zi, double& xscale, double& yscale, unsigned int width, unsigned int height, unsigned int frameno)
{
        zr = -0.743639266077433; 
        zi = +0.131824786875559; 
        double scale = 4. * std::pow(2, -std::min(frameno / 60., 53.)*0.7);
        xscale = scale / height; // TODO: Width??
        yscale = scale / height;
}

void setCoordinatesTime(double& zr, double& zi, double& xscale, double& yscale, unsigned int width, unsigned int height)
{
        zr = -0.743639266077433; 
        zi = +0.131824786875559; 
        double scale = 4. * std::pow(2, -std::min(GetTime(), 53.)*0.7);
        xscale = scale / height; // TODO: Width??
        yscale = scale / height;
}

std::tuple<float, float, float> getColor(unsigned int x, unsigned int y, float iter, bool process) {
    if(!process)
        return std::make_tuple(x, y, iter);
    else
    {
        static const unsigned char r[]{0x00,0x40,0x7E,0x11,0x16,0x38,0xFC,0xD0,0x5F,0xDC,0xFF,0x6B};
        static const unsigned char g[]{0x00,0x40,0x9F,0x90,0x68,0xCF,0xFF,0x99,0x00,0x37,0x8E,0x14};
        static const unsigned char b[]{0x00,0xE0,0xFF,0x9F,0x18,0x3F,0x00,0x24,0x09,0x0A,0xFE,0xBC};
        constexpr int k = 1, m = 0x3F;
		double d = ((((x&4)/4u + (x&2)*2u + (x&1)*16u) + (((x^y)&4)/2u + ((x^y)&2)*4u + ((x^y)&1)*32u))&m)/64.;
        auto lerp = [d,k](int a,int b,double p) -> unsigned { return int(a/k + (b/k-a/k) * p + d)*255/(255/k); };
        auto R = lerp(r[int(iter)%sizeof r], r[int(iter+1)%sizeof r], iter-int(iter));
        auto G = lerp(g[int(iter)%sizeof g], g[int(iter+1)%sizeof g], iter-int(iter));
        auto B = lerp(b[int(iter)%sizeof b], b[int(iter+1)%sizeof b], iter-int(iter));

		return std::make_tuple(R/255., G/255., B/255.);


    }
               
}


