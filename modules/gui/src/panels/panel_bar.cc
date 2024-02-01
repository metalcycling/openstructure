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


#include <ost/gui/gosty_app.hh>
#include <ost/gui/perspective.hh>
#include <ost/gui/widget.hh>
#include <ost/gui/widget_registry.hh>

#include "drop_box.hh"
#include "panel_bar.hh"
#include "panel_manager.hh"
#include "splitter_panel_bar.hh"
#include "tabbed_panel_bar.hh"

#include <QSettings>
#include <QVariant>
namespace ost {
namespace gui {

PanelBar::PanelBar(QWidget* parent) :
  Widget(NULL, parent),
  layout_(new QStackedLayout(this)),
  view_modes_( new QActionGroup(this)),
  view_mode_menu_(new QMenu("View Mode",this)),
  current_view_mode_(NULL),
  widget_states_(),
  drop_box_(new DropBox), // DropBox having no parent is done on purpose (otherwise it will not hide properly)
  show_action_(new QAction(this))
{
  connect(view_mode_menu_,SIGNAL(triggered(QAction*)),this,SLOT(ChangeViewMode(QAction*)));

  layout_->setMargin(0);
  layout_->setSpacing(0);
  this->setLayout(layout_);
  layout_->setCurrentIndex(0);

  this->setVisible(false);

  show_action_->setText("Show");
  show_action_->setCheckable(true);
  show_action_->setChecked(false);
  show_action_->setEnabled(false);
  connect(show_action_, SIGNAL(triggered(bool)), this,
          SLOT(ShowActionTrigger()));
  this->addAction(show_action_);
}

PanelBar::~PanelBar()
{
  delete drop_box_; // manually destroy drop_box_ widget, as it has no parent.
}

void PanelBar::AddWidget(Widget* widget, bool is_hidden)
{
  int index = this->GetIndex(widget);
  if (index == -1) {
    QString class_name = widget->GetUniqueID();
    WidgetState ws = WidgetState(widget, !is_hidden, class_name);
    widget_states_.append(ws);
    if(current_view_mode_){
      current_view_mode_->ShowWidget(widget, widget_states_.size() - 1, !is_hidden);
    }
    this->UpdateShow();
  }
}

void PanelBar::AddViewMode(PanelWidgetContainer * pwc)
{
  QAction* action = new QAction(pwc->GetName(),this);
  QVariant data = QVariant();
  data.setValue(pwc);
  action->setData(data);
  action->setCheckable(true);
  action->setChecked(false);
  view_modes_->addAction(action);
  view_mode_menu_->addAction(action);

  layout_->addWidget(pwc);
  if(!view_modes_->checkedAction()){
    action->setChecked(true);
  }
}

QMenu* PanelBar::GetViewModeMenu()
{
  return view_mode_menu_;
}

bool PanelBar::Save(const QString& prefix)
{
  QSettings settings;
  settings.beginGroup(prefix);
  settings.setValue("collapsed", this->isHidden());
  settings.setValue("type", layout_->currentIndex());
  settings.beginWriteArray("items");
  int index = 0;
  for (int i = 0; i < widget_states_.size(); ++i) {
    if (widget_states_.at(i).widget) {
      settings.setArrayIndex(i);
      settings.setValue("class_name", widget_states_.at(i).name);
      settings.setValue("show", widget_states_.at(i).visible);
      widget_states_.at(i).widget->Save(settings.group());
      ++index;
    }
  }
  settings.endArray();
  QList<QAction*> modes = view_modes_->actions();
  for(int i = 0; i< modes.size(); i++){
    PanelWidgetContainer* pwc = modes[i]->data().value<PanelWidgetContainer*>();
    if(pwc){
      pwc->Save(prefix);
    }
  }
  settings.endGroup(); // side_bar
  return true;
}

bool PanelBar::Restore(const QString& prefix)
{
  QSettings settings;
  settings.beginGroup(prefix);
  if (settings.contains("collapsed")) {
    this->setHidden(settings.value("collapsed").toBool());
    this->SetViewMode(settings.value("type").toInt());
    PanelPosition pos = NONE;
    if (settings.contains("position")){
      pos = PanelPosition(settings.value("position").toInt());
    }
    int size = settings.beginReadArray("items");
    for (int i = 0; i < size; ++i) {
      settings.setArrayIndex(i);
      QString class_name = settings.value("class_name").toString();
      bool show = settings.value("show").toBool();
      if(pos != NONE){
        GostyApp::Instance()->GetPerspective()->GetPanels()->AddWidgetByName(pos,class_name,!show);
      }
      else{
        Widget* widget = WidgetRegistry::Instance()->Create(class_name, NULL);
        if (widget == NULL) {
          continue;
        }
        this->AddWidget(widget, !show);
      }
    }
    settings.endArray();
  } else {
    this->SetViewMode(0);
    settings.endGroup();
    return false;
  }
  QList<QAction*> modes = view_modes_->actions();
  for(int i = 0; i< modes.size(); i++){
    PanelWidgetContainer* pwc = modes[i]->data().value<PanelWidgetContainer*>();
    if(pwc){
      pwc->Restore(prefix);
    }
  }

  for (int i = 0; i < widget_states_.size(); ++i) {
    if (widget_states_.at(i).widget) {
      QString path = prefix + "/items/" + QString::number(i+1);
      widget_states_.at(i).widget->Restore(path);
    }
  }
  settings.endGroup();
  return true;
}

void PanelBar::ToggleHide()
{
  show_action_->setChecked(!show_action_->isChecked());
  ShowActionTrigger();
}


QAction* PanelBar::GetShowAction(){
  return show_action_;
}

void PanelBar::ShowActionTrigger(){
  if (show_action_->isChecked()) {
      this->UpdateShow();
    } else {
      this->setHidden(true);
    }
}

void PanelBar::ToggleViewMode()
{
  int next_index = ((layout_->currentIndex() + 1) % layout_->count());
  this->SetViewMode(next_index);
}

void PanelBar::SetViewMode(int index)
{
  QList<QAction*> modes = view_modes_->actions();
  this->ChangeViewMode(modes[index]);
}

void PanelBar::InsertWidget(Widget * widget, int index, bool is_hidden)
{
  int oldindex = this->GetIndex(widget);
  if (oldindex >= 0) {
    widget_states_.move(oldindex, index);
    if(current_view_mode_){
      current_view_mode_->WidgetMoved(widget, index);
    }
  } else {
    QString class_name = widget->GetUniqueID();
    WidgetState ws = WidgetState(widget, !is_hidden, class_name);
    widget_states_.insert(index, ws);
    if(current_view_mode_){
      current_view_mode_->ShowWidget(widget_states_.at(index).widget, index, widget_states_.at(index).visible);
    }
  }
  this->UpdateShow();
}

void PanelBar::RemoveWidget(Widget *widget)
{
  if (widget) {
    int pos = -1;
    for (int i = 0; i < widget_states_.size(); i++) {
      if (widget_states_.at(i).widget == widget) {
        if(current_view_mode_){
          current_view_mode_->ShowWidget(widget_states_.at(i).widget, i, false);
        }
        pos = i;
        break;
      }
    }
    if (pos >= 0) {
      widget_states_.removeAt(pos);
    }
    this->UpdateShow();
  }
}

void PanelBar::ToggleVisibility(Widget* widget)
{
  int index = this->GetIndex(widget);
  if (index >= 0) {
    widget_states_.replace(index, WidgetState(widget,
        !widget_states_.at(index).visible, widget_states_.at(index).name));
    if(current_view_mode_){
      current_view_mode_->ShowWidget(widget, index, widget_states_.at(index).visible);
    }
  }
  this->UpdateShow();
}

void PanelBar::UpdateShow()
{
  bool show = false;
  for (int i = 0; i < widget_states_.size(); i++) {
    if (widget_states_.at(i).visible) {
      show = true;
      break;
    }
  }
  show_action_->setEnabled(show);
  show_action_->setChecked(show);
  this->setVisible(show);
}

int PanelBar::GetIndex(Widget* widget)
{
  for (int i = 0; i < widget_states_.size(); i++) {
    if (widget_states_.at(i).widget == widget) {
      return i;
    }
  }
  return -1;
}

void PanelBar::DragEnter(){
  if(this->isHidden()){
    if(current_view_mode_){
      current_view_mode_->StartDrag();
    }
    this->AddWidget(drop_box_,false);
    this->setVisible(true);
  }
}

void PanelBar::EndDrag(){
  this->RemoveWidget(drop_box_);
  if(current_view_mode_){
    current_view_mode_->EndDrag();
  }
}

void PanelBar::ChangeViewMode(QAction* action)
{
  if(action){
    PanelWidgetContainer* new_pwc = action->data().value<PanelWidgetContainer*>();
    if(new_pwc){
      if(current_view_mode_){
        for (int i = 0; i < widget_states_.size(); i++) {
          current_view_mode_->ShowWidget(widget_states_.at(i).widget, i, false);
        }
      }
      for (int i = 0; i < widget_states_.size(); i++) {
        new_pwc->ShowWidget(widget_states_.at(i).widget, i, widget_states_.at(i).visible);
      }
      layout_->setCurrentWidget(new_pwc);
      action->setChecked(true);
      current_view_mode_=new_pwc;
    }
  }
}

}}

