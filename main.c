#define SDL_MAIN_USE_CALLBACKS 1
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <math.h>

static SDL_Window *window;
static SDL_Renderer *renderer;
constexpr int weight = 1920, height = 1080;

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[]) {
  if (!SDL_Init(SDL_INIT_VIDEO))
    return SDL_APP_FAILURE;
  if (!SDL_CreateWindowAndRenderer("fft", weight, height, 0, &window,
                                   &renderer))
    return SDL_APP_FAILURE;
  return SDL_APP_CONTINUE;
}

SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event) {
  if (event->type == SDL_EVENT_QUIT)
    return SDL_APP_SUCCESS;
  return SDL_APP_CONTINUE;
}

float f(float x) {
  float f0 = 5e3, p0 = 0;
  return sinpi(f0 * x / 2 + p0);
}
SDL_AppResult SDL_AppIterate(void *appstate) {
  SDL_SetRenderDrawColorFloat(renderer, 0, 0, 0, 0);
  SDL_RenderClear(renderer);

  constexpr int n = 1024;
  float data[n];
  for (int i = 0; i < n; ++i)
    data[i] = f(i / 1e6);

  float h = height * .25;
  SDL_FPoint arr[n];
  for (int i = 0; i < n; ++i)
    arr[i].y = data[i] * -h + h, arr[i].x = i;

  SDL_SetRenderDrawColorFloat(renderer, 1, 1, 1, 1);
  SDL_RenderLines(renderer, arr, n);
  SDL_RenderPresent(renderer);
  return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result) {}
