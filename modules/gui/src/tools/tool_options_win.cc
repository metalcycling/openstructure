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
#include <ost/gui/tools/tool_options_win.hh>
#include <ost/gui/widget_registry.hh>
#include <ost/gui/tools/tool_options_widget.hh>
#include <ost/gui/tools/tool_manager.hh>

#include <QVBoxLayout>
namespace ost { namespace gui {


class ToolOptionsWinFactory : public WidgetFactory {
public:                                                                      
  ToolOptionsWinFactory(): WidgetFactory("ost::gui::ToolOptionsWin", 
                                         "Tool Options") {}

  Widget* Create(QWidget* parent)                                           
  {                                                                          
    return GostyApp::Instance()->GetToolOptionsWin();
  }                                                                          
};

OST_REGISTER_WIDGET(ToolOptionsWin, ToolOptionsWinFactory);

ToolOptionsWin::ToolOptionsWin(QWidget* parent):
  Widget(new QStackedWidget, parent)
{
  QVBoxLayout* l=new QVBoxLayout(this);
  l->setMargin(0);
  l->setSpacing(0);
  l->addWidget(this->GetInternalWidget());
  connect(&ToolManager::Instance(), SIGNAL(ToolAdded(const Tool*&)),
          this, SLOT(OnToolAdded(const Tool*&)));
  connect(&ToolManager::Instance(), SIGNAL(ActiveToolChanged(const Tool*&)),
          this, SLOT(OnActiveToolChange(const Tool*&)));
  const ToolList& tools=ToolManager::Instance().GetTools();
  for (ToolList::const_iterator i=tools.begin(), e=tools.end(); i!=e; ++i) {
    const Tool* tool=*i;
    this->OnToolAdded(const_cast<Tool*>(tool));
  }
  this->setMinimumSize(250,150);
}

void ToolOptionsWin::OnToolAdded(Tool* tool)
{
  ToolOptionsWidget* w=new ToolOptionsWidget(tool->GetToolOptions(), this);
  tool_to_widget_.insert(std::make_pair(tool, w));
  QStackedWidget* stacked;
  stacked=dynamic_cast<QStackedWidget*>(this->GetInternalWidget());
  stacked->addWidget(w);  
}

void ToolOptionsWin::OnActiveToolChange(Tool* tool)
{
  std::map<Tool*, ToolOptionsWidget*>::iterator i;
  i=tool_to_widget_.find(tool);
  assert(i!=tool_to_widget_.end());
  QStackedWidget* stacked;
  stacked=dynamic_cast<QStackedWidget*>(this->GetInternalWidget());
  stacked->setCurrentWidget(i->second);
}

bool ToolOptionsWin::Save(const QString& prefix)
{
  return true;
}

bool ToolOptionsWin::Restore(const QString& prefix)
{
  return true;
}

}}
