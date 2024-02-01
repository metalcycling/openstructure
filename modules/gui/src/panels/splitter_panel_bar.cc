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

#include "panel_widget_container.hh"
#include "splitter_panel_bar.hh"

#include <QSettings>

#include <QPalette>
namespace ost { namespace gui {

SplitterPanelBar::SplitterPanelBar(PanelBar* parent,Qt::Orientation orientation):
  PanelWidgetContainer(parent)
{
  layout_ = new QHBoxLayout(this);
  layout_->setMargin(0);
  layout_->setSpacing(0);  
  splitter_ = new ThinSplitter(orientation,this);
  sizes_ = QList<int>();
  layout_->addWidget(splitter_);

  connect(splitter_, SIGNAL(splitterMoved (int,int)), this,
          SLOT(SplitterMoved(int,int)));
}

bool SplitterPanelBar::Save(const QString& prefix)
{
  QSettings settings;
  settings.beginGroup(prefix + "/SplitterMode");
  settings.beginWriteArray("sizes");
  for (int i = 0; i < sizes_.size(); ++i) {
      settings.setArrayIndex(i);
      settings.setValue("length", sizes_.at(i));
  }
  settings.endArray();
  settings.endGroup();
  return true;
}



bool SplitterPanelBar::Restore(const QString& prefix)
{
QSettings settings;
  int size = settings.beginReadArray(prefix + "/SplitterMode/sizes");
  for (int i = 0; i < size; ++i) {
    settings.setArrayIndex(i);
    sizes_.append(settings.value("length").toInt());
    if(sizes_.size()>=splitter_->count()){
    splitter_->setSizes(sizes_);
    }
  }
  settings.endArray();
  return true;
}

void SplitterPanelBar::WidgetMoved(Widget* widget, int index){
  int pos=this->IndexOf(widget);
  if(pos >= 0){
    PanelBarWidgetHolder* wh = widget_holders_.at(pos);
    widget_holders_.removeAt(pos);
    wh->setVisible(false);
    wh->setParent(NULL);
  }
  PanelBarWidgetHolder* side_bar_widget_holder=new PanelBarWidgetHolder(widget,splitter_->orientation(),this->GetPanelBar());
  widget_holders_.insert(index,side_bar_widget_holder);
  splitter_->insertWidget(index, widget_holders_.at(index));
  if(pos>=0 && pos<sizes_.size() && index>=0 && index<sizes_.size()){
    sizes_.move(pos,index);
  }
  if(sizes_.size()>=splitter_->count()){
    splitter_->setSizes(sizes_);
  }
}


QString SplitterPanelBar::GetName()
{
  return "Splitter Mode";
}

void SplitterPanelBar::ShowWidget(Widget* widget, int pos, bool show){
  if(show)this->WidgetMoved(widget,pos);
  int index = this->IndexOf(widget);
  if(index>=0){
    widget_holders_.at(index)->GetWidget()->setVisible(show);
    widget_holders_.at(index)->setVisible(show);
    if(!show){
      widget_holders_.at(index)->setParent(NULL);
      widget_holders_.at(index)->close();
      widget_holders_.removeAt(index);
      sizes_ = splitter_->sizes();
    }
  }
}

int SplitterPanelBar::IndexOf(Widget* widget){
  for(int i = 0 ; i<widget_holders_.size();i++){
    if(widget_holders_.at(i)->GetWidget()==widget){
      return i;
    }
  }
  return -1;

}

void SplitterPanelBar::SplitterMoved(int pos, int index){
  sizes_ = splitter_->sizes();
}


SplitterPanelBar::~SplitterPanelBar(){
  //Do not destroy WidgetHolders, they might contain widgets which should not die x-(
  for(int i=0; i< widget_holders_.size(); i++){
    widget_holders_.at(i)->setParent(NULL);
  }
}

}}

