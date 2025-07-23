#include <SDL3/SDL_timer.h>
#include <complex.h>
#define SDL_MAIN_USE_CALLBACKS 1
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <fftw3.h>
#include <math.h>

static SDL_Window *window;
static SDL_Renderer *renderer;
constexpr int w = 2048, h = 1000, n = 1024;

static float *in;
static float complex *out;
static fftwf_plan p;
static SDL_FPoint arr[n];
SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[]) {
  if (!SDL_Init(SDL_INIT_VIDEO))
    return SDL_APP_FAILURE;
  if (!SDL_CreateWindowAndRenderer("fft", w, h, 0, &window, &renderer))
    return SDL_APP_FAILURE;
  in = fftwf_alloc_real(n);
  out = fftwf_alloc_complex(n);
  p = fftwf_plan_dft_r2c_1d(n, in, out, FFTW_MEASURE);
  for (int i = 0; i < n; ++i)
    arr[i].x = (float)i * w / n;
  return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event) {
  if (event->type == SDL_EVENT_QUIT)
    return SDL_APP_SUCCESS;
  return SDL_APP_CONTINUE;
}

static float f0 = 50e3, p0, f1 = 250e3, p1;
#define clk __builtin_readcyclecounter()
float f(float x) { return sinpi(f0 * x / 2 + p0) + sinpi(f1 * x / 2 + p1); }
SDL_AppResult SDL_AppIterate(void *appstate) {
  SDL_SetRenderDrawColorFloat(renderer, 0, 0, 0, 0);
  SDL_RenderClear(renderer);

  p0 = (clk >> 14 & 0xfffff) / 0xfffffp-1 - 1;
  p1 = (clk >> 18 & 0xfffff) / 0xfffffp-1 - .1;
  for (int i = 0; i < n; ++i)
    in[i] = f(i / 1e6);
  long t0 = clk;
  fftwf_execute(p);
  long t1 = clk;
  printf("%ld\n", t1 - t0);

  float y = h * .25;
  for (int i = 0; i < n; ++i)
    arr[i].y = in[i] * -y / 2 + y;
  SDL_SetRenderDrawColorFloat(renderer, 1, 1, 1, 1);
  SDL_RenderLines(renderer, arr, n);

  SDL_RenderDebugTextFormat(renderer, 0, 2 * y, "f0:%.0f", f0);
  SDL_RenderDebugTextFormat(renderer, w / 4., 2 * y, "f1:%.0f", f1);
  SDL_RenderDebugTextFormat(renderer, w / 2., 2 * y, "p0:%.5f", p0);
  SDL_RenderDebugTextFormat(renderer, w * .75, 2 * y, "p1:%.5f", p1);
  for (int i = 0; i < n; ++i)
    arr[i].y = crealf(out[i]) / n * -y + 3 * y;
  SDL_SetRenderDrawColorFloat(renderer, 0, 1, 1, 1);
  SDL_RenderLines(renderer, arr, n);
  for (int i = 0; i < n; ++i)
    arr[i].y = cimagf(out[i]) / n * -y + 3 * y;
  SDL_SetRenderDrawColorFloat(renderer, 1, 0, 1, 1);
  SDL_RenderLines(renderer, arr, n);

  SDL_RenderPresent(renderer);
  return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result) {}
