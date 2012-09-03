///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Dec 29 2008)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#ifndef __main_window_base__
#define __main_window_base__

#include <wx/intl.h>

#include <wx/panel.h>
#include <wx/gdicmn.h>
#include <wx/font.h>
#include <wx/colour.h>
#include <wx/settings.h>
#include <wx/string.h>
#include <wx/bitmap.h>
#include <wx/image.h>
#include <wx/icon.h>
#include <wx/bmpbuttn.h>
#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/sizer.h>
#include <wx/listctrl.h>
#include <wx/spinctrl.h>
#include <wx/slider.h>
#include <wx/statbox.h>
#include <wx/statline.h>
#include <wx/stattext.h>
#include <wx/choice.h>
#include <wx/splitter.h>
#include <wx/menu.h>
#include <wx/frame.h>

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
/// Class MainWindow_frame
///////////////////////////////////////////////////////////////////////////////
class MainWindow_frame : public wxFrame 
{
	private:
	
	protected:
		wxPanel* m_parent_panel1;
		wxPanel* m_panel1;
		wxBitmapButton* m_bpButton12;
		wxBitmapButton* m_bpButton132;
		wxBitmapButton* m_bpButton15;
		wxCheckBox* m_checkBox1;
		wxPanel* m_panel2;
		wxSplitterWindow* m_splitter3;
		wxPanel* m_panel3;
		wxListCtrl* m_listCtrl1;
		wxBitmapButton* m_bpButton13;
		wxBitmapButton* m_bpButton14;
		wxButton* m_button2;
		wxSpinCtrl* m_spinCtrl1;
		wxBitmapButton* m_bpButton1;
		wxBitmapButton* m_bpButton2;
		wxSlider* m_slider1;
		wxBitmapButton* m_bpButton3;
		wxBitmapButton* m_bpButton4;
		wxPanel* m_panel4;
		wxListCtrl* m_listCtrl11;
		wxBitmapButton* m_bpButton131;
		wxBitmapButton* m_bpButton141;
		wxBitmapButton* m_bpButton27;
		wxButton* m_button21;
		wxStaticLine* m_staticline1;
		wxStaticText* m_staticText2;
		wxStaticText* m_staticText3;
		wxChoice* m_choice1;
		wxBitmapButton* m_bpButton30;
		wxChoice* m_choice2;
		wxBitmapButton* m_bpButton31;
		wxMenuBar* m_menubar1;
		wxMenu* m_menu1;
		wxMenu* m_menu2;
		wxMenu* m_menu4;
		wxMenu* m_menu7;
		wxMenu* m_menu5;
		
		// Virtual event handlers, overide them in your derived class
		virtual void evt_reset_view( wxCommandEvent& event ){ event.Skip(); }
		virtual void evt_fit_to_window( wxCommandEvent& event ){ event.Skip(); }
		virtual void evt_background_color( wxCommandEvent& event ){ event.Skip(); }
		virtual void evt_perspective_toggle( wxCommandEvent& event ){ event.Skip(); }
		virtual void evt_system_add( wxCommandEvent& event ){ event.Skip(); }
		virtual void evt_system_delete( wxCommandEvent& event ){ event.Skip(); }
		virtual void evt_system_load( wxCommandEvent& event ){ event.Skip(); }
		
	
	public:
		MainWindow_frame( wxWindow* parent, wxWindowID id = wxID_ANY, const wxString& title = wxEmptyString, const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxSize( 707,560 ), long style = wxDEFAULT_FRAME_STYLE|wxTAB_TRAVERSAL );
		~MainWindow_frame();
		void m_splitter3OnIdle( wxIdleEvent& )
		{
		m_splitter3->SetSashPosition( 220 );
		m_splitter3->Disconnect( wxEVT_IDLE, wxIdleEventHandler( MainWindow_frame::m_splitter3OnIdle ), NULL, this );
		}
		
	
};

#endif //__main_window_base__
