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


#include <ost/gui/widget_registry.hh>

#include "bottom_bar.hh"
#include "side_bar.hh"
#include "panel_manager.hh"
#include "tabbed_panel_bar.hh"
#include "splitter_panel_bar.hh"

#include <QSettings>
#include <QMenu>
#include <QDebug>
namespace ost { namespace gui {

class ExternalWidgetFactory: public WidgetFactory {
private:
  Widget* widget_;
public:
  ExternalWidgetFactory(const QString& id, const QString& full_name, Widget* widget):
    WidgetFactory(id, full_name), widget_(widget) {
  }

  virtual Widget* Create(QWidget* parent) {
    return widget_;
  }
};


PanelManager::PanelManager(QWidget* widget):
 Widget(NULL, widget->parentWidget()),
     right_panel_splitter_(new ThinSplitter(Qt::Horizontal, this)),
     left_panel_splitter_(new ThinSplitter(Qt::Horizontal, this)),
     bottom_panel_splitter_(new ThinSplitter(Qt::Vertical, this)),
     pool_(new WidgetPool(this))
{

  PanelBar* bottom_bar = new BottomBar(this);
  PanelBar* right_bar = new SideBar(this);
  PanelBar* left_bar = new SideBar(this);

  right_panel_splitter_->addWidget(widget);
  right_panel_splitter_->addWidget(right_bar);

  bottom_panel_splitter_->addWidget(right_panel_splitter_);
  bottom_panel_splitter_->addWidget(bottom_bar);

  left_panel_splitter_->addWidget(left_bar);
  left_panel_splitter_->addWidget(bottom_panel_splitter_);
  //Set Default Width of left sidebar
  QList<int> left_panel_size;
  left_panel_size.append(210);
  left_panel_splitter_->setSizes(left_panel_size);

  // these lines make sure that the side bar remains constant when
  // resizing
  left_panel_splitter_->setStretchFactor(0, 0);
  left_panel_splitter_->setStretchFactor(1, 1);

  bottom_panel_splitter_->setStretchFactor(0, 0);
  bottom_panel_splitter_->setStretchFactor(1, 1);

  QHBoxLayout* layout = new QHBoxLayout(this);
  layout->setMargin(0);
  layout->setSpacing(0);
  layout->addWidget(left_panel_splitter_);
  this->setLayout(layout);

  panels_[LEFT_PANEL] = left_bar;
  panels_[RIGHT_PANEL] = right_bar;
  panels_[BOTTOM_PANEL] = bottom_bar;

  current_menu_ = NONE;
  add_menu_ = new QMenu("Add Widget", this);
  QAction* menu_action = add_menu_->menuAction();

  connect(menu_action,SIGNAL(hovered()),this,SLOT(UpdateAddMenu()));
  connect(add_menu_,SIGNAL(triggered(QAction*)),this,SLOT(AddWidgetFromMenu(QAction*)));
  this->acceptDrops();
}

void PanelManager::AddWidgetToPool(const QString& class_name, int limit){
  pool_->Add(class_name, limit);
}

void PanelManager::AddWidgetToPool(const QString& name, Widget* widget) {
  //widget->setParent(this);
  widget->SetUniqueID(name);
  ExternalWidgetFactory* ewf = new ExternalWidgetFactory(name, name, widget);
  WidgetRegistry::Instance()->RegisterWidgetFactory(ewf);
  pool_->Add(name, 1);
}

void PanelManager::AddWidget(PanelPosition pos, Widget* widget, bool is_hidden){
  if(panels_.contains(pos) && this->GetParentPanel(widget) == NONE){
    if(pool_->IsAvailable(widget->GetUniqueID())){
      panels_[pos]->AddWidget(widget,is_hidden);
      pool_->Take(widget->GetUniqueID());
    }
  }
}

void PanelManager::AddWidgetByName(PanelPosition pos, const QString& class_name,
                             bool is_hidden)
{
  QString n=class_name;
  n.replace(".", "::");
  if(panels_.contains(pos) && pool_->IsAvailable(n)) {
      WidgetRegistry* wf=WidgetRegistry::Instance();
      this->AddWidget(pos,wf->Create(n,this),is_hidden);
  }
}

bool PanelManager::Save(const QString& prefix)
{
  bool saveok = true;
  QSettings settings;
  settings.beginGroup(prefix);
  settings.beginGroup("left_panel");
  saveok = panels_[LEFT_PANEL]->Save(prefix+"/left_panel");
  settings.setValue("position",LEFT_PANEL);
  settings.beginWriteArray("sizes");
  QList<int> sizes = left_panel_splitter_->sizes();
  for (int i = 0; i < sizes.size(); ++i) {
      settings.setArrayIndex(i);
      settings.setValue("width", sizes.at(i));
  }
  settings.endArray();
  settings.endGroup();
  settings.beginGroup("bottom_panel");
  saveok = panels_[BOTTOM_PANEL]->Save(prefix+"/bottom_panel") && saveok;
  settings.setValue("position",BOTTOM_PANEL);
  settings.beginWriteArray("sizes");
  sizes = bottom_panel_splitter_->sizes();
  for (int i = 0; i < sizes.size(); ++i) {
      settings.setArrayIndex(i);
      settings.setValue("width", sizes.at(i));
  }
  settings.endArray();
  settings.endGroup();
  settings.beginGroup("right_panel");
  saveok = panels_[RIGHT_PANEL]->Save(prefix+"/right_panel") && saveok;
  settings.setValue("position",RIGHT_PANEL);
  settings.beginWriteArray("sizes");
  sizes = right_panel_splitter_->sizes();
  for (int i = 0; i < sizes.size(); ++i) {
      settings.setArrayIndex(i);
      settings.setValue("width", sizes.at(i));
  }
  settings.endArray();
  settings.endGroup();

  return saveok;
}

bool PanelManager::Restore(const QString& prefix)
{
  QList<int> sizes = QList<int>();
  QSettings settings;
  int size = settings.beginReadArray(prefix + "/left_panel/sizes");
  for (int i = 0; i < size; ++i) {
    settings.setArrayIndex(i);
    sizes.append(settings.value("width").toInt());
    if(sizes.size()>=left_panel_splitter_->count()){
      left_panel_splitter_->setSizes(sizes);
    }
  }
  settings.endArray();

  sizes.clear();
  size = settings.beginReadArray(prefix + "/bottom_panel/sizes");
  for (int i = 0; i < size; ++i) {
    settings.setArrayIndex(i);
    sizes.append(settings.value("width").toInt());
    if(sizes.size()>=bottom_panel_splitter_->count()){
      bottom_panel_splitter_->setSizes(sizes);
    }
  }
  settings.endArray();

  sizes.clear();
  size = settings.beginReadArray(prefix + "/right_panel/sizes");
  for (int i = 0; i < size; ++i) {
    settings.setArrayIndex(i);
    sizes.append(settings.value("width").toInt());
    if(sizes.size()>=right_panel_splitter_->count()){
      right_panel_splitter_->setSizes(sizes);
    }
  }
  settings.endArray();

  bool restoreok = panels_[LEFT_PANEL]->Restore(prefix+"/left_panel");
  restoreok = panels_[BOTTOM_PANEL]->Restore(prefix+"/bottom_panel") && restoreok;
  restoreok = panels_[RIGHT_PANEL]->Restore(prefix+"/right_panel") && restoreok;
  return restoreok;
}

void PanelManager::ToggleViewMode(PanelPosition pos){
  if(panels_.contains(pos)){
   panels_[pos]->ToggleViewMode();
  }
}

void PanelManager::ToggleHide(){
  QMapIterator<PanelPosition, PanelBar*> i(panels_);
  while (i.hasNext()) {
      i.next();
      this->ToggleHide(i.key());
  }
}

void PanelManager::ToggleHide(PanelPosition pos){
  if(panels_.contains(pos)){
   panels_[pos]->ToggleHide();
  }
}

void PanelManager::MoveWidget(Widget * widget, PanelPosition pos, int index){
  if(panels_.contains(pos)){
    PanelPosition current_position = this->GetParentPanel(widget);
    if(current_position == pos){
      panels_[pos]->InsertWidget(widget,index,false);
    }
    else{
      if(panels_.contains(current_position)){
        panels_[current_position]->RemoveWidget(widget);
        pool_->Give(widget->GetUniqueID());
      }
      if(pool_->IsAvailable(widget->GetUniqueID())){
        panels_[pos]->InsertWidget(widget,index,false);
        pool_->Take(widget->GetUniqueID());
      }
    }
  }
}

void PanelManager::MoveNextTo(Widget* target, Widget* widget){
  PanelPosition target_position =  this->GetParentPanel(target);
  if(panels_.contains(target_position)){
    this->MoveWidget(widget,target_position,panels_[target_position]->GetIndex(target));
  }
}

void PanelManager::ReplaceWidget(Widget* w1, Widget* w2){
  if(pool_->IsAvailable(w2->GetUniqueID())){
    PanelPosition pos_w1 = this->GetParentPanel(w1);
    PanelPosition pos_w2  = this->GetParentPanel(w2);
    if(panels_.contains(pos_w1) && pos_w2 == NONE){
      int index = panels_[pos_w1]->GetIndex(w1);
      panels_[pos_w1]->RemoveWidget(w1);
      pool_->Give(w1->GetUniqueID());
      panels_[pos_w1]->InsertWidget(w2,index,false);
      pool_->Take(w2->GetUniqueID());
    }
  }
}

void PanelManager::ReplaceWidget(Widget* w1, QString& class_name){
  if(pool_->IsAvailable(class_name)){
    PanelPosition pos = this->GetParentPanel(w1);
    if(panels_.contains(pos)){
      int index = panels_[pos]->GetIndex(w1);
      panels_[pos]->RemoveWidget(w1);
      pool_->Give(w1->GetUniqueID());
      WidgetRegistry* wf=WidgetRegistry::Instance();
      Widget* widget = wf->Create(class_name,this);
      panels_[pos]->InsertWidget(widget,index,false);
      pool_->Take(class_name);
    }
  }
}

void PanelManager::RemoveWidget(Widget *widget){
  PanelPosition pos = this->GetParentPanel(widget);
  if(panels_.contains(pos)){
    panels_[pos]->RemoveWidget(widget);
    pool_->Give(widget->GetUniqueID());
  }
}

QList<QString> PanelManager::GetAvailableWidgets()
{
  return pool_->GetAvailableWidgets();
}


QMenu* PanelManager::GetAvailableWidgetsMenu()
{
  QMenu* m=new QMenu;

  WidgetRegistry* wf=WidgetRegistry::Instance();
  QList<QString> available=this->GetAvailableWidgets();
  for (QList<QString>::iterator i=available.begin(),
       e=available.end(); i!=e; ++i) {
    QAction* action=m->addAction(wf->GetFullName(*i));
    action->setData(QVariant(*i));
  }
  return m;
}

PanelPosition PanelManager::GetParentPanel(Widget* widget){
  QMapIterator<PanelPosition, PanelBar*> i(panels_);
  while (i.hasNext()) {
      i.next();
      int index = i.value()->GetIndex(widget);
      if(index>=0){
        return i.key();
      }
  }
  return NONE;
}

void PanelManager::StartDrag()
{
  QMapIterator<PanelPosition, PanelBar*> i(panels_);
  while (i.hasNext()) {
    i.next();
    i.value()->DragEnter();
  }
}

void PanelManager::EndDrag()
{
  QMapIterator<PanelPosition, PanelBar*> i(panels_);
  while (i.hasNext()) {
    i.next();
    i.value()->EndDrag();
  }
}

QMenu* PanelManager::GetMenu()
{
  QMenu* panel_menu = new QMenu(this);
  panel_menu->setTitle("Panels");

  QMenu* left_panel = panel_menu->addMenu("Left Panel");
  QAction* menu_action = left_panel->menuAction();
  menu_action->setData(QVariant(LEFT_PANEL));
  connect(menu_action,SIGNAL(hovered()),this,SLOT(UpdateCurrentMenu()));
  left_panel->addAction(panels_[LEFT_PANEL]->GetShowAction());
  left_panel->addMenu(add_menu_);
  left_panel->addMenu(panels_[LEFT_PANEL]->GetViewModeMenu());

  QMenu* bottom_panel = panel_menu->addMenu("Bottom Panel");
  menu_action = bottom_panel->menuAction();
  menu_action->setData(QVariant(BOTTOM_PANEL));
  connect(menu_action,SIGNAL(hovered()),this,SLOT(UpdateCurrentMenu()));
  bottom_panel->addAction(panels_[BOTTOM_PANEL]->GetShowAction());
  bottom_panel->addMenu(add_menu_);
  bottom_panel->addMenu(panels_[BOTTOM_PANEL]->GetViewModeMenu());

  QMenu* right_panel = panel_menu->addMenu("Right Panel");
  menu_action = right_panel->menuAction();
  menu_action->setData(QVariant(RIGHT_PANEL));
  connect(menu_action,SIGNAL(hovered()),this,SLOT(UpdateCurrentMenu()));
  right_panel->addAction(panels_[RIGHT_PANEL]->GetShowAction());
  right_panel->addMenu(add_menu_);
  right_panel->addMenu(panels_[RIGHT_PANEL]->GetViewModeMenu());
  return panel_menu;
}

void PanelManager::UpdateAddMenu(){
  add_menu_->clear();
  WidgetRegistry* wf=WidgetRegistry::Instance();
  QList<QString> available=this->GetAvailableWidgets();
  for (QList<QString>::iterator i=available.begin(),
       e=available.end(); i!=e; ++i) {
    QAction* action=add_menu_->addAction(wf->GetFullName(*i));
    action->setData(QVariant(*i));
  }
}

void PanelManager::AddWidgetFromMenu(QAction* action){
  if(action){
    QString new_class_name=action->data().toString();
    if(panels_.contains(current_menu_)){
      this->AddWidgetByName(current_menu_,new_class_name);
    }
  }
}

void PanelManager::UpdateCurrentMenu(){
  QAction* action = qobject_cast<QAction*>(this->sender());
  if(action){
    current_menu_ = PanelPosition(action->data().toInt());
  }
}

}}

