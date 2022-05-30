/**
 * @file main.cpp  
 *
 * @brief GTI320 Labo 3 - Simulation d'un système masse-ressort
 *
 * Nom: William Lebel
 * Email : william.lebel.1@ens.etsmtl.ca
 *
 */

#include <nanogui/opengl.h>
#include <nanogui/glutil.h>
#include <nanogui/screen.h>
#include <nanogui/formhelper.h>
#include <nanogui/layout.h>

#include "ParticleSimApplication.h"

// Includes for the GLTexture class.
#include <ctime>

#if defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
#if defined(_WIN32)
#  pragma warning(push)
#  pragma warning(disable: 4457 4456 4005 4312)
#endif

#if defined(_WIN32)
#  pragma warning(pop)
#endif
#if defined(_WIN32)
#  if defined(APIENTRY)
#    undef APIENTRY
#  endif
#  include <windows.h>
#endif


int main(int argc, char** argv) 
{
  nanogui::init();
  nanogui::ref<ParticleSimApplication> app = new ParticleSimApplication();
  app->drawAll();
  app->setVisible(true);
  nanogui::mainloop(4);
  nanogui::shutdown();

  return 0;
}
