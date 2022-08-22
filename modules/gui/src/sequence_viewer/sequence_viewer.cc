//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------

/*
  Author: Stefan Scheuber
 */
#include <ost/gui/gosty_app.hh>
#include <boost/pointer_cast.hpp>


#include <ost/platform.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/entity_view.hh>

#include <ost/seq/sequence_handle.hh>

#include <ost/gfx/entity.hh>
#include <ost/gfx/scene.hh>
#include <ost/gfx/gfx_node_visitor.hh>

#include <ost/gui/widget_registry.hh>

#include "sequence_model.hh"
#include "sequence_viewer.hh"
#include "sequence_search_bar.hh"
#include "sequence_table_view.hh"

#include <QAbstractItemView>
#include <QApplication>
#include <QClipboard>
#include <QDir>
#include <QHeaderView>
#include <QMenu>
#include <QPushButton>
#include <QShortcut>
#include <QVBoxLayout>
#include <QVarLengthArray>

namespace ost { namespace gui {

class SequenceViewerFactory: public WidgetFactory {
public:
  SequenceViewerFactory() :
    WidgetFactory("ost::gui::SequenceViewer", "Sequence Viewer") {
  }

  virtual Widget* Create(QWidget* parent) {
    return GostyApp::Instance()->GetSequenceViewer();
  }
};

OST_REGISTER_WIDGET(SequenceViewer, SequenceViewerFactory);

struct NodeParentPair {
  NodeParentPair(gfx::GfxNodeP n, gfx::GfxNodeP p): node(n), parent(p) {}
  
  gfx::GfxNodeP node;
  gfx::GfxNodeP parent;
};

struct GetNodesVisitor: public gfx::GfxNodeVisitor {

  GetNodesVisitor(): nodes_() {}
  virtual void VisitObject(gfx::GfxObj* o, const Stack& st) {
    nodes_.push_back(o->shared_from_this());    
  }
  gfx::GfxNodeVector nodes_;
  gfx::GfxNodeVector GetNodes() {return nodes_; }
};

SequenceViewer::SequenceViewer(bool stand_alone, bool observe_scene,
                               const QString& title,
                               QWidget* parent): Widget(NULL,parent)
{
  observe_scene_=observe_scene;
  model_ = new SequenceModel(this);
  this->setWindowTitle(title);
  QVBoxLayout* layout = new QVBoxLayout(this);
  layout->setMargin(0);
  layout->setSpacing(0);
  this->setLayout(layout);

  toolbar_ = new QToolBar(this);
  this->InitActions();


  if(stand_alone){
    this->InitMenuBar();
  }
  this->InitSearchBar();
  this->InitView();
  if (observe_scene_) {
    gfx::Scene::Instance().AttachObserver(this);
    gfx::GfxNodeP root_node = gfx::Scene::Instance().GetRootNode();
    GetNodesVisitor gnv;
    gfx::Scene::Instance().Apply(gnv);
    gfx::GfxNodeVector list=gnv.GetNodes();
    for(unsigned int i=0; i<list.size();i++){
      this->NodeAdded(list[i]);
    }    
  }
  connect(seq_table_view_, SIGNAL(AlignmentChanged()),
          this, SIGNAL(AlignmentChanged()));
}

void SequenceViewer::InitMenuBar()
{
  toolbar_->setToolButtonStyle(Qt::ToolButtonIconOnly);
  toolbar_->setIconSize(QSize(16,16));
  toolbar_->addActions(action_list_);
  layout()->addWidget(toolbar_);
}

void SequenceViewer::InitSearchBar()
{
  seq_search_bar_ = new SeqSearchBar(this);
  seq_search_bar_->hide();
  layout()->addWidget(seq_search_bar_);
  connect(seq_search_bar_, SIGNAL(Changed(const QString&, bool, const QString&)), 
         this, SLOT(OnSearchBarUpdate(const QString&, bool, const QString&)));
}

void SequenceViewer::InitView()
{
  seq_table_view_ = new SequenceTableView(model_);
  seq_table_view_->setFrameShape(QFrame::NoFrame);
  seq_table_view_->setTextElideMode(Qt::ElideRight);  
  layout()->addWidget(seq_table_view_);

  connect(model_,SIGNAL(columnsInserted(const QModelIndex&, int, int)),seq_table_view_,SLOT(columnCountChanged(const QModelIndex&, int, int)));
  connect(model_,SIGNAL(rowsInserted(const QModelIndex&, int, int)),seq_table_view_,SLOT(rowCountChanged(const QModelIndex&, int, int)));

  seq_table_view_->horizontalHeader()->setMinimumSectionSize(2);
  seq_table_view_->verticalHeader()->setMinimumSectionSize(2);
  seq_table_view_->setSelectionMode(QAbstractItemView::ExtendedSelection);
  connect(seq_table_view_->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this, SLOT(SelectionModelChanged(const QItemSelection&, const QItemSelection&)));
  connect(seq_table_view_,SIGNAL(doubleClicked(const QModelIndex&)),model_,SLOT(DoubleClicked(const QModelIndex&)));
#if !(defined(__APPLE__) && (QT_VERSION>=0x040600))  
  connect(seq_table_view_->GetStaticColumn(),SIGNAL(doubleClicked(const QModelIndex&)),this,SLOT(DoubleClicked(const QModelIndex&)));
  connect(seq_table_view_->GetStaticRow(),SIGNAL(doubleClicked(const QModelIndex&)),this,SLOT(DoubleClicked(const QModelIndex&)));
#endif
  connect(seq_table_view_,SIGNAL(CopyEvent(QKeyEvent*)),this,SLOT(CopyEvent(QKeyEvent*)));
  connect(seq_table_view_,SIGNAL(MouseWheelEvent(QWheelEvent*)),this,SLOT(MouseWheelEvent(QWheelEvent*)));
}

void SequenceViewer::InitActions()
{
  QDir icon_path(GetSharedDataPath().c_str());
  icon_path.cd("gui");
  icon_path.cd("icons");

  QAction* find_action = new QAction(toolbar_);
  find_action->setText("Find Dialog");
  find_action->setShortcut(QKeySequence(tr("Ctrl+F")));
  find_action->setCheckable(true);
  find_action->setToolTip("Display Find-Dialog (Ctrl+F)");
  find_action->setIcon(QIcon(icon_path.absolutePath()+QDir::separator()+QString("find_icon.png")));
  action_list_.append(find_action);
  connect(find_action, SIGNAL(triggered(bool)), this, SLOT(FindInSequence()));

  display_mode_actions_ = new QActionGroup(toolbar_);
  QAction* menu_action = new QAction(toolbar_);
  menu_action->setText("Menubar");
  menu_action->setShortcut(QKeySequence(tr("Ctrl+M")));
  menu_action->setToolTip("Display Options (Ctrl+M)");
  menu_action->setIcon(QIcon(icon_path.absolutePath()+QDir::separator()+QString("menubar_icon.png")));
  action_list_.append(menu_action);
  connect(menu_action, SIGNAL(triggered(bool)), this, SLOT(DisplayMenu()));
}

void SequenceViewer::AddEntity(const gfx::EntityP& entity)
{
  if (!observe_scene_) {
    return;
  }
  seq_table_view_->SetSingleAlignment(false);
  model_->InsertGfxEntity(entity);
  this->FitToContents();
  this->UpdateSearchBar();
}

void SequenceViewer::RemoveEntity(const gfx::EntityP& entity)
{
  if (!observe_scene_) {
    return;
  }
  model_->RemoveGfxEntity(entity);
}

void SequenceViewer::NodeAdded(const gfx::GfxNodeP& n)
{
  if (!observe_scene_) {
    return;
  }
  if (gfx::EntityP o=boost::dynamic_pointer_cast<gfx::Entity>(n)) {
    this->AddEntity(o);
  }
}

void SequenceViewer::NodeRemoved(const gfx::GfxNodeP& node)
{
  if (!observe_scene_) {
    return;
  }
  if (gfx::EntityP o=boost::dynamic_pointer_cast<gfx::Entity>(node)) {
    this->RemoveEntity(o);
  }
}

void SequenceViewer::AddAlignment(const seq::AlignmentHandle& alignment)
{
  if(alignment.GetCount()>0 && alignment.GetLength()>0){
    if (model_->rowCount()>1) {
      seq_table_view_->SetSingleAlignment(false);
    } else {
      seq_table_view_->SetSingleAlignment(true);
    }
    model_->InsertAlignment(alignment);
    this->FitToContents();
  }
}

void SequenceViewer::SetAlignment(const seq::AlignmentHandle& alignment)
{
  model_->Clear();
  this->AddAlignment(alignment);
}

void SequenceViewer::RemoveAlignment(const seq::AlignmentHandle& alignment)
{
  model_->RemoveAlignment(alignment);
}

void SequenceViewer::UpdateSearchBar()
{
  QStringList sequence_names_;
  for(int i = 1; i< model_->rowCount(); i++){
    QString name = model_->data(model_->index(i,0),Qt::DisplayRole).toString();
    sequence_names_.append(name);
  }
  seq_search_bar_->UpdateItems(sequence_names_);
}

void SequenceViewer::SelectionModelChanged(const QItemSelection& sel, const QItemSelection& desel)
{
  gfx::Scene::Instance().DetachObserver(this);
  model_->SelectionChanged(sel, desel);
  gfx::Scene::Instance().AttachObserver(this);
}

void SequenceViewer::SelectionChanged(const gfx::GfxObjP& o,
                                      const mol::EntityView& view)
{
  disconnect(seq_table_view_->selectionModel(), 
             SIGNAL(selectionChanged(const QItemSelection&, 
                                     const QItemSelection&)), 
             this, SLOT(SelectionModelChanged(const QItemSelection&, const QItemSelection&)));
  gfx::EntityP entity=boost::dynamic_pointer_cast<gfx::Entity>(o);
  if(entity){
    const QModelIndexList& list = model_->GetModelIndexes(entity, view);
    this->SelectList(list);
  }
  connect(seq_table_view_->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this, SLOT(SelectionModelChanged(const QItemSelection&, const QItemSelection&)));
}

void SequenceViewer::DoubleClicked(const QModelIndex& index)
{
  disconnect(seq_table_view_->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this, SLOT(SelectionModelChanged(const QItemSelection&, const QItemSelection&)));
  model_->DoubleClicked(index);
  connect(seq_table_view_->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this, SLOT(SelectionModelChanged(const QItemSelection&, const QItemSelection&)));
}

void SequenceViewer::FitToContents()
{
  seq_table_view_->resizeColumnsToContents();
  seq_table_view_->resizeRowsToContents();
  seq_table_view_->setColumnWidth(0, 150);  
}

void SequenceViewer::MouseWheelEvent(QWheelEvent* event)
{
  int delta = event->angleDelta().y();
  if (event->orientation() == Qt::Vertical) {
    if(delta>0){
      model_->ZoomIn();
      this->FitToContents();
      seq_table_view_->viewport()->update();      
    }
    else if(delta<0){
      model_->ZoomOut();
      this->FitToContents();
      seq_table_view_->viewport()->update();      
    }
  }
  event->accept();
}

void SequenceViewer::CopyEvent(QKeyEvent* event)
{
  QItemSelectionModel* model = seq_table_view_->selectionModel();
  const QModelIndexList& list = model->selectedIndexes();
  if(! list.empty()){
    QString clipboard_string;
    QSet<int> rows;
    int min_col=model_->columnCount();
    int max_col=0;
    for(int i = 0; i < list.size(); i++){
      if(list[i].column()>max_col){
        max_col = list[i].column();
      }
      if(list[i].column()<min_col){
        min_col = list[i].column();
      }
      rows.insert(list[i].row());
    }

    bool first_row = true;
    for(int i = 1; i < model_->rowCount(); i++){
      if(rows.contains(i)){
        if(!first_row){
          clipboard_string.append("\n");
        }
        for(int j=min_col; j<=max_col; j++){
          const QModelIndex& index = model_->index(i,j);
          if(model->isSelected(index)){
            clipboard_string.append(model_->data(index,Qt::DisplayRole).toString());
          }
          else{
            clipboard_string.append('-');
          }
        }
        first_row = false;
      }
    }
    QApplication::clipboard()->setText(clipboard_string);
  }
  event->accept();
}

void SequenceViewer::FindInSequence()
{
  if(seq_search_bar_->isHidden()){
    seq_search_bar_->show();
  }
  else{
    seq_search_bar_->hide();
  }
}

void SequenceViewer::OnSearchBarUpdate(const QString& subject,
                                           bool search_in_all, const QString& name)
{
  seq_table_view_->selectionModel()->clear();
  if(search_in_all){
    const QModelIndexList& list = model_->GetModelIndexes(subject);
    this->SelectList(list);
  }
  else{
    const QModelIndexList& list = model_->GetModelIndexes(subject,name);
    this->SelectList(list);
  }
}

void SequenceViewer::SelectList(const QModelIndexList& list)
{
  QItemSelectionModel* model = seq_table_view_->selectionModel();
  QSet<int> rows_visited;
  for(int i = 0; i<list.size(); i++){
    int row =list[i].row();
    if(!rows_visited.contains(row)){
      model->select(list[i],QItemSelectionModel::Rows|QItemSelectionModel::Deselect);
      rows_visited.insert(row);
    }
  }
  if (! list.empty()) {
    int last_row = 0;
    int last_col = 0;
    QModelIndex topleft_idx;
    QItemSelection *selection = new QItemSelection();
    int i = 1;
    topleft_idx = list[0];
    last_row = list[0].row();
    last_col = list[0].column();
    for (i = 1; i < list.size(); i++) {
      // store block on discontinued row or unequal column
      if (((last_col + 1) != list[i].column()) || (last_row != list[i].row())) {
        selection->select(topleft_idx, list[i-1]);
        model->select(*selection, QItemSelectionModel::Select);
        topleft_idx = list[i];
      }
      last_row = list[i].row();
      last_col = list[i].column();
    }
    // store last block
    selection->select(topleft_idx, list[i-1]);
    model->select(*selection, QItemSelectionModel::Select);
  }
}

const QStringList& SequenceViewer::GetDisplayModes()
{
  return model_->GetDisplayModes();
}
const QStringList& SequenceViewer::GetDisplayModes(const seq::AlignmentHandle& alignment)
{
  return model_->GetDisplayModes(alignment);
}
const QStringList& SequenceViewer::GetDisplayModes(const gfx::EntityP& entity)
{
  return model_->GetDisplayModes(entity);
}

const QString& SequenceViewer::GetCurrentDisplayMode()
{
  return model_->GetCurrentDisplayMode();
}
const QString& SequenceViewer::GetCurrentDisplayMode(const seq::AlignmentHandle& alignment)
{
  return model_->GetCurrentDisplayMode(alignment);
}
const QString& SequenceViewer::GetCurrentDisplayMode(const gfx::EntityP& entity)
{
  return model_->GetCurrentDisplayMode(entity);
}

void SequenceViewer::ChangeDisplayMode(const QString& string)
{
  model_->SetDisplayMode(string);
  seq_table_view_->viewport()->update();
}

void SequenceViewer::ChangeDisplayMode(const seq::AlignmentHandle& alignment, const QString& string)
{
  model_->SetDisplayMode(alignment, string);
  seq_table_view_->viewport()->update();
}

void SequenceViewer::ChangeDisplayMode(const gfx::EntityP& entity, const QString& string)
{
  model_->SetDisplayMode(entity, string);
  seq_table_view_->viewport()->update();
}

ActionList SequenceViewer::GetActions()
{
  return action_list_;
}

void SequenceViewer::DisplayMenu()
{
  QMenu* menu = new QMenu();
  QList<QAction*> actions = display_mode_actions_->actions();
  for(int i=0;i<actions.size();i++){
    display_mode_actions_->removeAction(actions[i]);
  }
  const QStringList& display_modes = this->GetDisplayModes();
  for(int i=0; i<display_modes.size(); i++){
    QString ident(display_modes[i]);
    QAction* action = new QAction(ident,menu);
    action->setCheckable(true);
    connect(action,SIGNAL(triggered(bool)),this,SLOT(ChangeDisplayMode()));
    display_mode_actions_->addAction(action);
    if(display_modes[i] == this->GetCurrentDisplayMode() ){
      action->setChecked(true);
    }
    menu->addAction(action);
  }
  menu->exec(QCursor::pos());
}

void SequenceViewer::ChangeDisplayMode()
{
  QAction* action = display_mode_actions_->checkedAction();
  if(action){
    this->ChangeDisplayMode(action->text());
  }
}

SequenceViewer::~SequenceViewer(){
  gfx::Scene::Instance().DetachObserver(this);
}

}}
