#include <SDL3/SDL_timer.h>
#include <complex.h>
#include <stdio.h>
#define SDL_MAIN_USE_CALLBACKS 1
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <math.h>

static SDL_Window *window;
static SDL_Renderer *renderer;
constexpr int w = 2048, h = 1000, n = 64, logn = __builtin_ctz(n);

static SDL_FPoint arr[w > n ? w : n];

typedef float f32;
typedef complex float c32;
static int rev[n];
static f32 in[n];
static c32 out[n], rot[n];
void fft_init() {
  for (int i = 0; i < n; ++i)
    rev[i] = __builtin_bitreverse32(i) >> (32 - logn),
    rot[i] = cexp(-2i * M_PI * i / n);
}
void fft(f32 const *restrict in, c32 *restrict out) {
  for (int i = 0; i < n; ++i)
    out[rev[i]] = in[i];
  for (int i = 0; i < logn; ++i) {
    int m = 1 << i;
    for (int j = 0; j < n; j += m << 1)
      for (int k = 0; k < m; ++k) {
        c32 t = rot[k << (logn - i - 1)] * out[j + k + m];
        c32 u = out[j + k];
        out[j + k] = u + t;
        out[j + k + m] = u - t;
      }
  }
}

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[]) {
  if (!SDL_Init(SDL_INIT_VIDEO))
    return SDL_APP_FAILURE;
  if (!SDL_CreateWindowAndRenderer("fft", w, h, 0, &window, &renderer))
    return SDL_APP_FAILURE;
  fft_init();
  return SDL_APP_CONTINUE;
}

static f32 f0 = 50e3, p0, f1 = 250e3, p1;
f32 tri(f32 x) { return fabs(2 - fmodf(x * 2, 4)) - 1; }
f32 f(f32 x) { return sinpi(f0 * 2 * x + p0) + tri(f1 * 2 * x + p1); }
SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event) {
  switch (event->type) {
  case SDL_EVENT_QUIT:
    return SDL_APP_SUCCESS;
  case SDL_EVENT_MOUSE_MOTION: {
    auto tmp = event->motion;
    f32 x = tmp.x / w, y = 1 - tmp.y / h;
    if (tmp.state != 1)
      goto end;
    if (x < .5) {
      if (x < .25)
        f0 = y * 5e5;
      else
        f1 = y * 5e5;
    } else {
      if (x < .75)
        p0 = y * 2 - 1;
      else
        p1 = y * 2 - 1;
    }
  }
  }
end:
  return SDL_APP_CONTINUE;
}

#define clk __builtin_readcyclecounter()
SDL_AppResult SDL_AppIterate(void *appstate) {
  SDL_SetRenderDrawColorFloat(renderer, 0, 0, 0, 1);
  SDL_RenderClear(renderer);

  for (int i = 0; i < n; ++i)
    in[i] = f(i / 1e6);

  static long st, cnt;
  ++cnt;
  long t0 = clk;
  fft(in, out);
  long t1 = clk;
  st += t1 - t0;
  printf("%f\n", st / (f32)cnt);

  f32 y = h * .25;
  for (int i = 0; i < w; ++i)
    arr[i].x = i;
  for (int i = 0; i < w; ++i)
    arr[i].y = f(i / 1e6 * n / w) * -y / 2 + y;
  SDL_SetRenderDrawColorFloat(renderer, .5, .5, .5, 1);
  SDL_RenderLines(renderer, arr, w);

  for (int i = 0; i < n; ++i)
    arr[i].x = (f32)i * w / n;
  for (int i = 0; i < n; ++i)
    arr[i].y = in[i] * -y / 2 + y;
  SDL_SetRenderDrawColorFloat(renderer, 1, 1, 1, .5);
  SDL_RenderLines(renderer, arr, n);

  SDL_RenderDebugTextFormat(renderer, 0, 2 * y, "f0:%.0f", f0);
  SDL_RenderDebugTextFormat(renderer, w / 4., 2 * y, "f1:%.0f", f1);
  SDL_RenderDebugTextFormat(renderer, w / 2., 2 * y, "p0:%.5f", p0);
  SDL_RenderDebugTextFormat(renderer, w * .75, 2 * y, "p1:%.5f", p1);

  for (int i = 0; i < n; ++i)
    arr[i].y = cabsf(out[i]) / n * 2 * -y + h;
  SDL_SetRenderDrawColorFloat(renderer, 0, 1, 0, .5);
  SDL_RenderLines(renderer, arr, n);
  for (int i = 0; i < n; ++i)
    arr[i].y = cargf(out[i]) / M_PI * -y + 3 * y;
  SDL_SetRenderDrawColorFloat(renderer, 1, 0, 0, .5);
  SDL_RenderLines(renderer, arr, n);

  SDL_RenderPresent(renderer);
  return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result) {}
