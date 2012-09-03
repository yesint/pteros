#ifndef __WXWIDGETSAPP_H
#define __WXWIDGETSAPP_H

#include <wx/wx.h>
#include "pteros/pteros_visualizer.h"

class App : public wxApp{
public:
    App();
    virtual ~App();
    virtual bool OnInit();
    pteros::System sys;
    pteros::Selection sel;

    // Pointer to the Visuzlizer_wx object
    pteros::Visualizer_wx* vis;
};

DECLARE_APP(App)

#endif





