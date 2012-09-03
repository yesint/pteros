#include "app.h"
#include "main_window.h"

IMPLEMENT_APP(App)

App::App(){
}

App::~App()
{
}

bool App::OnInit()
{
    // Load image handlers by hand, otherwise this dumb thing doesn't work
    wxInitAllImageHandlers();

    main_window* main_w = new main_window( (wxWindow*)NULL );
    main_w ->Show();
    SetTopWindow( main_w );

    // Set pointer to visualizer
    vis = (pteros::Visualizer_wx*) main_w->FindWindow(wxT("GLCanvas"));

    sys.load("2lao.pdb");
    sel = pteros::Selection(sys,"all");
    vis->show(sel);
    vis->center_all();

    return true;
}

