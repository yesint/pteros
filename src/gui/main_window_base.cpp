///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Dec 29 2008)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#include <wx/wxprec.h>

#include "main_window_base.h"

///////////////////////////////////////////////////////////////////////////

MainWindow_frame::MainWindow_frame( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxFrame( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );
	
	wxBoxSizer* bSizer1;
	bSizer1 = new wxBoxSizer( wxHORIZONTAL );
	
	m_parent_panel1 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER|wxTAB_TRAVERSAL );
	wxBoxSizer* bSizer10;
	bSizer10 = new wxBoxSizer( wxVERTICAL );
	
	m_panel1 = new wxPanel( m_parent_panel1, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSUNKEN_BORDER|wxTAB_TRAVERSAL );
	bSizer10->Add( m_panel1, 1, wxEXPAND|wxALL, 1 );
	
	wxBoxSizer* bSizer12;
	bSizer12 = new wxBoxSizer( wxHORIZONTAL );
	
	m_bpButton12 = new wxBitmapButton( m_parent_panel1, wxID_ANY, wxBitmap( wxT("icons/stock_home.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton12->SetToolTip( _("Reser view") );
	
	m_bpButton12->SetToolTip( _("Reser view") );
	
	bSizer12->Add( m_bpButton12, 0, wxALL, 5 );
	
	m_bpButton132 = new wxBitmapButton( m_parent_panel1, wxID_ANY, wxBitmap( wxT("icons/stock_fullscreen.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton132->SetToolTip( _("Fit selected to window") );
	
	m_bpButton132->SetToolTip( _("Fit selected to window") );
	
	bSizer12->Add( m_bpButton132, 0, wxALL, 5 );
	
	m_bpButton15 = new wxBitmapButton( m_parent_panel1, wxID_ANY, wxBitmap( wxT("icons/applications-graphics.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton15->SetToolTip( _("Background color") );
	
	m_bpButton15->SetToolTip( _("Background color") );
	
	bSizer12->Add( m_bpButton15, 0, wxALL, 5 );
	
	m_checkBox1 = new wxCheckBox( m_parent_panel1, wxID_ANY, _("Perspective"), wxDefaultPosition, wxDefaultSize, 0 );
	
	bSizer12->Add( m_checkBox1, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5 );
	
	bSizer10->Add( bSizer12, 0, wxEXPAND, 5 );
	
	m_parent_panel1->SetSizer( bSizer10 );
	m_parent_panel1->Layout();
	bSizer10->Fit( m_parent_panel1 );
	bSizer1->Add( m_parent_panel1, 1, wxEXPAND|wxTOP|wxBOTTOM|wxLEFT, 5 );
	
	m_panel2 = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	wxBoxSizer* bSizer2;
	bSizer2 = new wxBoxSizer( wxVERTICAL );
	
	m_splitter3 = new wxSplitterWindow( m_panel2, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxSP_3D );
	m_splitter3->Connect( wxEVT_IDLE, wxIdleEventHandler( MainWindow_frame::m_splitter3OnIdle ), NULL, this );
	m_panel3 = new wxPanel( m_splitter3, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxHSCROLL|wxSUNKEN_BORDER|wxTAB_TRAVERSAL|wxVSCROLL );
	wxStaticBoxSizer* sbSizer2;
	sbSizer2 = new wxStaticBoxSizer( new wxStaticBox( m_panel3, wxID_ANY, _("Systems") ), wxVERTICAL );
	
	sbSizer2->SetMinSize( wxSize( -1,100 ) ); 
	m_listCtrl1 = new wxListCtrl( m_panel3, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLC_EDIT_LABELS|wxLC_HRULES|wxLC_REPORT );
	sbSizer2->Add( m_listCtrl1, 1, wxEXPAND|wxTOP|wxBOTTOM, 5 );
	
	wxBoxSizer* bSizer4;
	bSizer4 = new wxBoxSizer( wxVERTICAL );
	
	wxBoxSizer* bSizer6;
	bSizer6 = new wxBoxSizer( wxHORIZONTAL );
	
	m_bpButton13 = new wxBitmapButton( m_panel3, wxID_ANY, wxBitmap( wxT("icons/list-add.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton13->SetToolTip( _("New system") );
	
	m_bpButton13->SetToolTip( _("New system") );
	
	bSizer6->Add( m_bpButton13, 0, wxALIGN_CENTER_VERTICAL|wxBOTTOM|wxRIGHT, 5 );
	
	m_bpButton14 = new wxBitmapButton( m_panel3, wxID_ANY, wxBitmap( wxT("icons/list-remove.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton14->SetToolTip( _("Delete system") );
	
	m_bpButton14->SetToolTip( _("Delete system") );
	
	bSizer6->Add( m_bpButton14, 0, wxALIGN_CENTER_VERTICAL|wxBOTTOM|wxRIGHT|wxLEFT, 5 );
	
	m_button2 = new wxButton( m_panel3, wxID_ANY, _("Load..."), wxDefaultPosition, wxDefaultSize, 0 );
	m_button2->SetToolTip( _("Load data to selected system") );
	
	bSizer6->Add( m_button2, 0, wxALIGN_CENTER_VERTICAL|wxBOTTOM|wxRIGHT|wxLEFT, 5 );
	
	bSizer4->Add( bSizer6, 0, wxEXPAND, 5 );
	
	wxBoxSizer* bSizer5;
	bSizer5 = new wxBoxSizer( wxHORIZONTAL );
	
	m_spinCtrl1 = new wxSpinCtrl( m_panel3, wxID_ANY, wxEmptyString, wxDefaultPosition, wxSize( 50,-1 ), wxSP_ARROW_KEYS, 0, 10, 0 );
	bSizer5->Add( m_spinCtrl1, 0, wxALIGN_CENTER_VERTICAL|wxTOP|wxBOTTOM|wxRIGHT, 5 );
	
	m_bpButton1 = new wxBitmapButton( m_panel3, wxID_ANY, wxBitmap( wxT("icons/go-first.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	bSizer5->Add( m_bpButton1, 0, wxALIGN_CENTER_VERTICAL|wxTOP|wxBOTTOM, 5 );
	
	m_bpButton2 = new wxBitmapButton( m_panel3, wxID_ANY, wxBitmap( wxT("icons/go-previous.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	bSizer5->Add( m_bpButton2, 0, wxALIGN_CENTER_VERTICAL|wxTOP|wxBOTTOM, 5 );
	
	m_slider1 = new wxSlider( m_panel3, wxID_ANY, 50, 0, 100, wxDefaultPosition, wxSize( 100,-1 ), wxSL_HORIZONTAL );
	bSizer5->Add( m_slider1, 1, wxALIGN_CENTER_VERTICAL|wxTOP|wxBOTTOM, 5 );
	
	m_bpButton3 = new wxBitmapButton( m_panel3, wxID_ANY, wxBitmap( wxT("icons/go-next.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	bSizer5->Add( m_bpButton3, 0, wxALIGN_CENTER_VERTICAL|wxTOP|wxBOTTOM, 5 );
	
	m_bpButton4 = new wxBitmapButton( m_panel3, wxID_ANY, wxBitmap( wxT("icons/go-last.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	bSizer5->Add( m_bpButton4, 0, wxALIGN_CENTER_VERTICAL|wxTOP|wxBOTTOM, 5 );
	
	bSizer4->Add( bSizer5, 0, wxEXPAND, 5 );
	
	sbSizer2->Add( bSizer4, 0, wxEXPAND, 5 );
	
	m_panel3->SetSizer( sbSizer2 );
	m_panel3->Layout();
	sbSizer2->Fit( m_panel3 );
	m_panel4 = new wxPanel( m_splitter3, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxHSCROLL|wxSUNKEN_BORDER|wxTAB_TRAVERSAL|wxVSCROLL );
	wxStaticBoxSizer* sbSizer21;
	sbSizer21 = new wxStaticBoxSizer( new wxStaticBox( m_panel4, wxID_ANY, _("Selections") ), wxVERTICAL );
	
	sbSizer21->SetMinSize( wxSize( -1,100 ) ); 
	m_listCtrl11 = new wxListCtrl( m_panel4, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLC_HRULES|wxLC_REPORT );
	sbSizer21->Add( m_listCtrl11, 1, wxEXPAND|wxTOP, 5 );
	
	wxBoxSizer* bSizer41;
	bSizer41 = new wxBoxSizer( wxVERTICAL );
	
	wxBoxSizer* bSizer61;
	bSizer61 = new wxBoxSizer( wxHORIZONTAL );
	
	m_bpButton131 = new wxBitmapButton( m_panel4, wxID_ANY, wxBitmap( wxT("icons/list-add.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton131->SetToolTip( _("New selection") );
	
	m_bpButton131->SetToolTip( _("New selection") );
	
	bSizer61->Add( m_bpButton131, 0, wxTOP|wxBOTTOM|wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	
	m_bpButton141 = new wxBitmapButton( m_panel4, wxID_ANY, wxBitmap( wxT("icons/list-remove.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton141->SetToolTip( _("Delete selection") );
	
	m_bpButton141->SetToolTip( _("Delete selection") );
	
	bSizer61->Add( m_bpButton141, 0, wxALIGN_CENTER_VERTICAL|wxTOP|wxBOTTOM|wxRIGHT, 5 );
	
	m_bpButton27 = new wxBitmapButton( m_panel4, wxID_ANY, wxBitmap( wxT("icons/gtk-edit.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton27->SetToolTip( _("Edit selection") );
	
	m_bpButton27->SetToolTip( _("Edit selection") );
	
	bSizer61->Add( m_bpButton27, 0, wxALL, 5 );
	
	m_button21 = new wxButton( m_panel4, wxID_ANY, _("Save..."), wxDefaultPosition, wxDefaultSize, 0 );
	m_button21->SetToolTip( _("Save selection to file") );
	
	bSizer61->Add( m_button21, 0, wxALIGN_CENTER_VERTICAL|wxTOP|wxBOTTOM, 5 );
	
	bSizer41->Add( bSizer61, 0, wxEXPAND, 5 );
	
	m_staticline1 = new wxStaticLine( m_panel4, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL );
	bSizer41->Add( m_staticline1, 0, wxEXPAND|wxRIGHT|wxLEFT, 5 );
	
	wxFlexGridSizer* fgSizer1;
	fgSizer1 = new wxFlexGridSizer( 2, 2, 0, 0 );
	fgSizer1->SetFlexibleDirection( wxBOTH );
	fgSizer1->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	
	m_staticText2 = new wxStaticText( m_panel4, wxID_ANY, _("Representation"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText2->Wrap( -1 );
	fgSizer1->Add( m_staticText2, 0, wxTOP|wxRIGHT|wxLEFT, 5 );
	
	m_staticText3 = new wxStaticText( m_panel4, wxID_ANY, _("Coloring"), wxDefaultPosition, wxDefaultSize, 0 );
	m_staticText3->Wrap( -1 );
	fgSizer1->Add( m_staticText3, 0, wxTOP|wxRIGHT|wxLEFT, 5 );
	
	wxBoxSizer* bSizer21;
	bSizer21 = new wxBoxSizer( wxHORIZONTAL );
	
	wxString m_choice1Choices[] = { _("None"), _("CPK"), _("VDW"), _("Lines"), _("Bonds") };
	int m_choice1NChoices = sizeof( m_choice1Choices ) / sizeof( wxString );
	m_choice1 = new wxChoice( m_panel4, wxID_ANY, wxDefaultPosition, wxDefaultSize, m_choice1NChoices, m_choice1Choices, 0 );
	m_choice1->SetSelection( 0 );
	bSizer21->Add( m_choice1, 0, wxTOP|wxBOTTOM|wxLEFT|wxALIGN_CENTER_VERTICAL, 5 );
	
	m_bpButton30 = new wxBitmapButton( m_panel4, wxID_ANY, wxBitmap( wxT("icons/gtk-edit.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton30->SetToolTip( _("Edit representation") );
	
	m_bpButton30->SetToolTip( _("Edit representation") );
	
	bSizer21->Add( m_bpButton30, 0, wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	
	fgSizer1->Add( bSizer21, 1, wxEXPAND, 5 );
	
	wxBoxSizer* bSizer23;
	bSizer23 = new wxBoxSizer( wxHORIZONTAL );
	
	wxString m_choice2Choices[] = { _("By element"), _("By beta"), _("By occupancy"), _("Single color") };
	int m_choice2NChoices = sizeof( m_choice2Choices ) / sizeof( wxString );
	m_choice2 = new wxChoice( m_panel4, wxID_ANY, wxDefaultPosition, wxDefaultSize, m_choice2NChoices, m_choice2Choices, 0 );
	m_choice2->SetSelection( 0 );
	bSizer23->Add( m_choice2, 0, wxALIGN_CENTER_VERTICAL|wxLEFT, 5 );
	
	m_bpButton31 = new wxBitmapButton( m_panel4, wxID_ANY, wxBitmap( wxT("icons/gtk-edit.png"), wxBITMAP_TYPE_ANY ), wxDefaultPosition, wxDefaultSize, wxBU_AUTODRAW );
	m_bpButton31->SetToolTip( _("Edit coloring") );
	
	m_bpButton31->SetToolTip( _("Edit coloring") );
	
	bSizer23->Add( m_bpButton31, 0, wxRIGHT|wxALIGN_CENTER_VERTICAL, 5 );
	
	fgSizer1->Add( bSizer23, 1, wxEXPAND, 5 );
	
	bSizer41->Add( fgSizer1, 1, wxEXPAND, 5 );
	
	sbSizer21->Add( bSizer41, 0, wxEXPAND, 5 );
	
	m_panel4->SetSizer( sbSizer21 );
	m_panel4->Layout();
	sbSizer21->Fit( m_panel4 );
	m_splitter3->SplitHorizontally( m_panel3, m_panel4, 220 );
	bSizer2->Add( m_splitter3, 1, wxEXPAND, 5 );
	
	m_panel2->SetSizer( bSizer2 );
	m_panel2->Layout();
	bSizer2->Fit( m_panel2 );
	bSizer1->Add( m_panel2, 0, wxEXPAND | wxALL, 5 );
	
	this->SetSizer( bSizer1 );
	this->Layout();
	m_menubar1 = new wxMenuBar( 0 );
	m_menu1 = new wxMenu();
	wxMenuItem* m_menuItem1;
	m_menuItem1 = new wxMenuItem( m_menu1, wxID_ANY, wxString( _("New project...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu1->Append( m_menuItem1 );
	
	wxMenuItem* m_menuItem2;
	m_menuItem2 = new wxMenuItem( m_menu1, wxID_ANY, wxString( _("Load project...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu1->Append( m_menuItem2 );
	
	wxMenuItem* m_menuItem3;
	m_menuItem3 = new wxMenuItem( m_menu1, wxID_ANY, wxString( _("Save project...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu1->Append( m_menuItem3 );
	
	m_menubar1->Append( m_menu1, _("Project") );
	
	m_menu2 = new wxMenu();
	wxMenuItem* m_menuItem4;
	m_menuItem4 = new wxMenuItem( m_menu2, wxID_ANY, wxString( _("New") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu2->Append( m_menuItem4 );
	
	wxMenuItem* m_menuItem6;
	m_menuItem6 = new wxMenuItem( m_menu2, wxID_ANY, wxString( _("Delete") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu2->Append( m_menuItem6 );
	
	wxMenuItem* m_menuItem5;
	m_menuItem5 = new wxMenuItem( m_menu2, wxID_ANY, wxString( _("Load data...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu2->Append( m_menuItem5 );
	
	wxMenuItem* m_menuItem8;
	m_menuItem8 = new wxMenuItem( m_menu2, wxID_ANY, wxString( _("Delete frames...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu2->Append( m_menuItem8 );
	
	wxMenuItem* m_menuItem9;
	m_menuItem9 = new wxMenuItem( m_menu2, wxID_ANY, wxString( _("Delete all selections") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu2->Append( m_menuItem9 );
	
	m_menubar1->Append( m_menu2, _("System") );
	
	m_menu4 = new wxMenu();
	wxMenuItem* m_menuItem10;
	m_menuItem10 = new wxMenuItem( m_menu4, wxID_ANY, wxString( _("New...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu4->Append( m_menuItem10 );
	
	wxMenuItem* m_menuItem11;
	m_menuItem11 = new wxMenuItem( m_menu4, wxID_ANY, wxString( _("Delete") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu4->Append( m_menuItem11 );
	
	wxMenuItem* m_menuItem12;
	m_menuItem12 = new wxMenuItem( m_menu4, wxID_ANY, wxString( _("Edit...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu4->Append( m_menuItem12 );
	
	m_menu4->AppendSeparator();
	
	wxMenuItem* m_menuItem25;
	m_menuItem25 = new wxMenuItem( m_menu4, wxID_ANY, wxString( _("Representation...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu4->Append( m_menuItem25 );
	
	wxMenuItem* m_menuItem26;
	m_menuItem26 = new wxMenuItem( m_menu4, wxID_ANY, wxString( _("Coloring...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu4->Append( m_menuItem26 );
	
	m_menubar1->Append( m_menu4, _("Selection") );
	
	m_menu7 = new wxMenu();
	wxMenuItem* m_menuItem14;
	m_menuItem14 = new wxMenuItem( m_menu7, wxID_ANY, wxString( _("Translate...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu7->Append( m_menuItem14 );
	
	wxMenuItem* m_menuItem15;
	m_menuItem15 = new wxMenuItem( m_menu7, wxID_ANY, wxString( _("Rotate...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu7->Append( m_menuItem15 );
	
	wxMenuItem* m_menuItem16;
	m_menuItem16 = new wxMenuItem( m_menu7, wxID_ANY, wxString( _("Fit") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu7->Append( m_menuItem16 );
	
	wxMenuItem* m_menuItem17;
	m_menuItem17 = new wxMenuItem( m_menu7, wxID_ANY, wxString( _("RMSD") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu7->Append( m_menuItem17 );
	
	wxMenuItem* m_menuItem18;
	m_menuItem18 = new wxMenuItem( m_menu7, wxID_ANY, wxString( _("Fit trajectory") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu7->Append( m_menuItem18 );
	
	wxMenuItem* m_menuItem24;
	m_menuItem24 = new wxMenuItem( m_menu7, wxID_ANY, wxString( _("RMSD trajectory") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu7->Append( m_menuItem24 );
	
	m_menubar1->Append( m_menu7, _("Actions") );
	
	m_menu5 = new wxMenu();
	wxMenuItem* m_menuItem20;
	m_menuItem20 = new wxMenuItem( m_menu5, wxID_ANY, wxString( _("Perspective projection") ) , wxEmptyString, wxITEM_CHECK );
	m_menu5->Append( m_menuItem20 );
	
	wxMenuItem* m_menuItem21;
	m_menuItem21 = new wxMenuItem( m_menu5, wxID_ANY, wxString( _("Reset view") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu5->Append( m_menuItem21 );
	
	wxMenuItem* m_menuItem22;
	m_menuItem22 = new wxMenuItem( m_menu5, wxID_ANY, wxString( _("Fit selected to window") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu5->Append( m_menuItem22 );
	
	wxMenuItem* m_menuItem23;
	m_menuItem23 = new wxMenuItem( m_menu5, wxID_ANY, wxString( _("Background color...") ) , wxEmptyString, wxITEM_NORMAL );
	m_menu5->Append( m_menuItem23 );
	
	m_menubar1->Append( m_menu5, _("View") );
	
	this->SetMenuBar( m_menubar1 );
	
	
	// Connect Events
	m_bpButton12->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_reset_view ), NULL, this );
	m_bpButton132->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_fit_to_window ), NULL, this );
	m_bpButton15->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_background_color ), NULL, this );
	m_checkBox1->Connect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_perspective_toggle ), NULL, this );
	m_bpButton13->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_system_add ), NULL, this );
	m_bpButton14->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_system_delete ), NULL, this );
	m_button2->Connect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_system_load ), NULL, this );
	this->Connect( m_menuItem4->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_system_add ) );
	this->Connect( m_menuItem6->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_system_delete ) );
	this->Connect( m_menuItem5->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_system_load ) );
	this->Connect( m_menuItem20->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_perspective_toggle ) );
	this->Connect( m_menuItem21->GetId(), wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_reset_view ) );
}

MainWindow_frame::~MainWindow_frame()
{
	// Disconnect Events
	m_bpButton12->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_reset_view ), NULL, this );
	m_bpButton132->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_fit_to_window ), NULL, this );
	m_bpButton15->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_background_color ), NULL, this );
	m_checkBox1->Disconnect( wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_perspective_toggle ), NULL, this );
	m_bpButton13->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_system_add ), NULL, this );
	m_bpButton14->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_system_delete ), NULL, this );
	m_button2->Disconnect( wxEVT_COMMAND_BUTTON_CLICKED, wxCommandEventHandler( MainWindow_frame::evt_system_load ), NULL, this );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_system_add ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_system_delete ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_system_load ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_perspective_toggle ) );
	this->Disconnect( wxID_ANY, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler( MainWindow_frame::evt_reset_view ) );
}
