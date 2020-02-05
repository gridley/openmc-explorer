#include <SDL2/SDL.h>
#include <openmc/position.h>
#include <openmc/error.h>
#include <openmc/math_functions.h>
#include <array>
#include <string>
#include <cmath>
#include <cstdlib>
#include <iostream>
/* TODO
  - Use timers to steady the scrolling rate
  - Make TraceableGeometry into an abstract base class
    that defines how Monte Carlo codes should integrate
    with this program
*/

// In actuality, use openmc::Position
struct Ray
{
  openmc::Position r;
  openmc::Direction d;
};

std::array<unsigned char, 3> hsv2rgb(float h, float s, float l) {
  std::array<unsigned char,3> result;
  std::array<float,3> resultf;
  float c = (1. - abs(2.*l - 1.)) * s;
  int hprime = h / 60.;
  float x = c * (1.-abs(fmod(hprime,2)-1));
  switch(hprime) {
    case 0:
      resultf[0]=c;
      resultf[1]=x;
      resultf[2]=0;
      break;
    case 1:
      resultf[0]=x;
      resultf[1]=c;
      resultf[2]=0;
      break;
    case 2:
      resultf[0]=0;
      resultf[1]=c;
      resultf[2]=x;
      break;
    case 3:
      resultf[0]=0;
      resultf[1]=x;
      resultf[2]=c;
      break;
    case 4:
      resultf[0]=x;
      resultf[1]=0;
      resultf[2]=c;
      break;
    case 5:
      resultf[0]=c;
      resultf[1]=0;
      resultf[2]=x;
      break;
  }
  float m = l - c/2;
  result[0] = static_cast<unsigned char>(255*(resultf[0]+m));
  result[1] = static_cast<unsigned char>(255*(resultf[1]+m));
  result[2] = static_cast<unsigned char>(255*(resultf[2]+m));
  return result;
}

// Abstract this eventually
class TraceableGeometry
{
  double cube_rad;
  public:
    TraceableGeometry();
    double traceToNextSurface(Ray& ray);
    int getMaterialID();
};

TraceableGeometry::TraceableGeometry() : cube_rad(1.0) {}

double
TraceableGeometry::traceToNextSurface(Ray& ray)
{
  // Check intersections with each surface
  std::array<double, 6> intersect_distances;
}

// Simple struct that holds settings for the program, and has default values
struct ExplorerSettings
{
  bool running;
  double delta_angle; // degrees to rotate camera on keypress
  double delta_scoot; // distance to scoot cammera (cm)

  // How many lines to draw between frames. This dynamically adjusts to keep
  // a desirable framerate
  unsigned n_lines_to_draw = 1;
  
  unsigned window_size_x;
  unsigned window_size_y;
  unsigned n_threads();

  ExplorerSettings();
};
ExplorerSettings::ExplorerSettings() : running(true), delta_angle(5.0), delta_scoot(5.),
  n_lines_to_draw(5), window_size_x(1920), window_size_y(1080), n_threads(1) {}

class WindowManager
{
  public:
    int window_width, window_height;
    SDL_Window* window;
    WindowManager(ExplorerSettings& sett);
    ~WindowManager();
};
WindowManager::WindowManager(ExplorerSettings& sett) : window_width(sett.window_size_x),
  window_height(sett.window_size_y),
  window(nullptr)
{
  SDL_Init(SDL_INIT_VIDEO);
  window = SDL_CreateWindow("OpenMC Explorer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
      window_width, window_height, SDL_WINDOW_SHOWN);

  if (window == nullptr)
    openmc::fatal_error("Failed to create window instance.\n");
}
WindowManager::~WindowManager()
{
  SDL_DestroyWindow(window);
  SDL_Quit();
}

class Camera
{
  openmc::Position pos;
  openmc::Direction dir;

  // whether to update the rendering
  bool viewChanged_;

  // Amount of roll to the camera (radians)
  double roll;

  // Field of view. Initially set to values close to the human eye
  double fov_x;
  double fov_y;

  // Needs a window manager object to know how many rays to create
  WindowManager& window;

  public:

    // Initialize at (1,1,1) pointing at origin
    Camera(WindowManager& win);

    // Phi and mu are given in degrees here:
    Camera(double x, double y, double z, double phi, double mu, WindowManager& win);

    void rotatePhi(double angle_d);
    void rotateMu(double angle_d);
    void increment_roll(double angle_d);
    void increment_fov_horiz(double angle_d);
    void increment_fov_vert(double angle_d);

    // Scoot in a Cartesian direction
    template<unsigned Dir>
    void scoot(double distance) {
      pos[Dir] += distance;
      viewChanged_ = true;
    }

    // Check whether view has been changed (resets internal state after checking)
    bool viewChanged();
};

Camera::Camera(double x, double y, double z, double phi, double mu,
    WindowManager& win) : dir(x, y, z),
  fov_x(135.0*M_PI/180.),
  fov_y(75.0*M_PI/180.),
  roll(0.0),
  window(win)
{
  double phi_rad = phi * M_PI/180.0;
  double mu_rad = mu * M_PI/180.0;
  dir.x = cos(phi_rad) * sin(mu_rad);
  dir.y = sin(phi_rad) * sin(mu_rad);
  dir.z = cos(mu_rad);
}
Camera::Camera(WindowManager& win) : Camera(1.0, 1.0, 1.0, 225.0, 135.0, win) {}

void Camera::rotatePhi(double angle_d) {
  double angle_r = angle_d * M_PI/180.0;
  dir = openmc::rotate_angle(dir, 0.0, &angle_r, nullptr);
  viewChanged_ = true;
}

void Camera::rotateMu(double angle_d) {
  double phi_zero = 0.0;
  double angle_r = angle_d * M_PI/180.;
  dir = openmc::rotate_angle(dir, angle_r, &phi_zero, nullptr);
  viewChanged_ = true;
}
bool Camera::viewChanged() {
  bool result = viewChanged_;
  viewChanged_ = false;
  return result;
}
void Camera::increment_roll(double angle_d) {
  roll += angle_d * M_PI/180.;
  viewChanged_ = true;
}
void Camera::increment_fov_horiz(double angle_d) {
  fov_x += angle_d * M_PI/180.;
  viewChanged_ = true;
}
void Camera::increment_fov_vert(double angle_d) {
  fov_y += angle_d * M_PI/180.;
  viewChanged_ = true;
}

void handle_events(Camera& cam, ExplorerSettings& sett) {
  SDL_Event event;
  // Process keyboard input
  while (SDL_PollEvent(&event)) {
    switch (event.type) {
      case SDL_QUIT:
        sett.running=false;
        break;
      case SDL_KEYDOWN:
        switch (event.key.keysym.sym) {

          // Camera rotation
          case SDLK_UP:
            cam.rotateMu(sett.delta_angle);
            break;
          case SDLK_DOWN:
            cam.rotateMu(-sett.delta_angle);
            break;
          case SDLK_LEFT:
            cam.rotatePhi(-sett.delta_angle);
            break;
          case SDLK_RIGHT:
            cam.rotatePhi(sett.delta_angle);
            break;

          // Camera movement
          case SDLK_w:
            cam.scoot<0>(sett.delta_scoot);
            break;
          case SDLK_s:
            cam.scoot<0>(-sett.delta_scoot);
            break;
          case SDLK_a:
            cam.scoot<1>(-sett.delta_scoot);
            break;
          case SDLK_d:
            cam.scoot<1>(sett.delta_scoot);
            break;
          case SDLK_LSHIFT:
            cam.scoot<2>(-sett.delta_scoot);
            break;
          case SDLK_SPACE:
            cam.scoot<2>(sett.delta_scoot);
            break;

          // Camera roll
          case SDLK_e:
            cam.increment_roll(sett.delta_angle);
            break;
          case SDLK_q:
            cam.increment_roll(-sett.delta_angle);
            break;

          // Field of view adjustment
          case SDLK_r:
            cam.increment_fov_horiz(sett.delta_angle);
            break;
          case SDLK_f:
            cam.increment_fov_horiz(-sett.delta_angle);
            break;
          case SDLK_t:
            cam.increment_fov_vert(sett.delta_angle);
            break;
          case SDLK_g:
            cam.increment_fov_vert(-sett.delta_angle);
            break;
        }
        break;
    }
  }
}

// Simple class to do RAII for SDL renderers
class Renderer
{
  SDL_Renderer* renderer_;
  public:
    Renderer(WindowManager& win);
    ~Renderer();
    void present();
};
Renderer::Renderer(WindowManager& win) : renderer_(nullptr) {
  renderer_ = SDL_CreateRenderer(win.window, -1, SDL_RENDERER_ACCELERATED || SDL_HINT_RENDER_VSYNC);
  if (renderer_ == nullptr)
  {
    openmc::fatal_error("unable to make renderer instance: \n"
        + std::string(SDL_GetError()));
  }
  // Initialize black screen
  SDL_RenderClear(renderer_);
  SDL_SetRenderDrawColor(renderer_, 0, 0, 0, 0);
  SDL_RenderClear(renderer_);
}
void Renderer::present() { SDL_RenderPresent(renderer_); }
Renderer::~Renderer() { SDL_DestroyRenderer(renderer_); }

int main(int argc, char* argv)
{
  ExplorerSettings settings;
  WindowManager window(settings);
  Renderer renderer(window);
  Camera cam(window);

  // Set up openmc
  int err;
  err = openmc_init(argc, argv, nullptr);
  if (err == -1) { return 0; }
  else if (err) { openmc::fatal_error(openmc_err_msg); }

  while (settings.running) {
    handle_events(cam, settings);

    if (cam.viewChanged() || cam.isRendering())
    {
      if (not cam.isRendering()) cam.beginRendering();
      cam.renderlines(renderer, sett)
    }

    renderer.present();
  }
}
