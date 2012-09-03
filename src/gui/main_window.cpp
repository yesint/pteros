#include "main_window.h"

using namespace pteros;

main_window::main_window( wxWindow* parent )
    :MainWindow_frame(parent)
{
    // Create sizer
    wxBoxSizer* sizer = new wxBoxSizer(wxHORIZONTAL);
    // Create visualizer
    Visualizer_wx* canvas = new Visualizer_wx(m_panel1);
    // Add visualizer to sizer
    sizer->Add(canvas, 1, wxEXPAND);
    // pack sizer to enclosing panel
    m_panel1->SetSizer(sizer);
    m_panel1->SetAutoLayout(true);

    // Set pointer to visualizer
    vis = canvas;
}

void main_window::evt_perspective_toggle( wxCommandEvent& event ){
	vis->camera.toggle_perspective();
	vis->Refresh();
}


void main_window::evt_fit_to_window( wxCommandEvent& event ){
    cout << "Centering in window" << endl;
    vis->reset_view();
    vis->Refresh();
}

void main_window::evt_system_add( wxCommandEvent& event )
{
	// TODO: Implement evt_system_add
}

void main_window::evt_system_delete( wxCommandEvent& event )
{
	// TODO: Implement evt_system_delete
}

void main_window::evt_system_load( wxCommandEvent& event )
{
	// TODO: Implement evt_system_load
}
