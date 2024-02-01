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
  Author: Marco Biasini
 */

#include "python_shell_widget.hh"
#include "python_shell.hh"
#include <ost/gui/gosty_app.hh>
#include <ost/gui/widget_registry.hh>


#include <QVBoxLayout>
namespace ost { namespace gui {

class PythonShellFactory: public WidgetFactory {
public:
	PythonShellFactory() :
    WidgetFactory("ost::gui::PythonShell", "Python Shell") {
  }

  virtual Widget* Create(QWidget* parent) {
    return GostyApp::Instance()->GetPyShell();
  }
};

OST_REGISTER_WIDGET(PythonShell, PythonShellFactory);


PythonShell::PythonShell(QWidget* parent):
  Widget(new PythonShellWidget(NULL), parent)
{
  QVBoxLayout* l=new QVBoxLayout(this);
  l->setSpacing(0);
  l->setMargin(0);
  l->addWidget(this->PyShell());
}

bool PythonShell::Restore(const QString& prefix)
{
  return true;
}

bool PythonShell::Save(const QString& prefix)
{
  return true;  
}


PythonShellWidget* PythonShell::PyShell()
{
  return dynamic_cast<PythonShellWidget*>(this->GetInternalWidget());
}

}}
