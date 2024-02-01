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
#include "tool_options_widget.hh"

#include <QGridLayout>
#include <QLabel>
#include <QPalette>
#include <QComboBox>
#include <QLineEdit>
#include <QIntValidator>
#include <QDoubleValidator>
#include <QSpacerItem>
#include <QPushButton>

namespace ost { namespace gui {

ToolOptionsWidget::ToolOptionsWidget(ToolOptions* options, 
                                     QWidget* parent):
  QWidget(parent), options_(options)
{
  this->Populate();
}

QWidget* ToolOptionsWidget::MakeIntWidget(ToolOptionInt* opts)
{
  QLineEdit* line=new QLineEdit(this);
  QIntValidator* vali=new QIntValidator(line);
  vali->setRange(opts->GetLowerLimit(), opts->GetUpperLimit());
  line->setValidator(vali);  
  QString t;
  t.setNum(opts->GetValue());
  line->setText(QString::number(opts->GetValue()));
  connect(line, SIGNAL(textChanged(const QString&)), 
          new IntOptBinder(line, opts), SLOT(TextChanged(const QString&)));
  return line;
}

QWidget* ToolOptionsWidget::MakeEnumWidget(ToolOptionEnum* opts)
{
  QComboBox* combo=new QComboBox(this);
  for (ToolOptionEnum::ConstIterator i=opts->Begin(), 
       e=opts->End(); i!=e; ++i) {
      combo->addItem(QString(i->text.c_str()), QVariant(i->tag));
  }
  connect(combo, SIGNAL(currentIndexChanged(int)), 
          new EnumOptBinder(combo, opts), SLOT(IndexChanged(int)));
  return combo;
} 
 
QWidget* ToolOptionsWidget::MakeFloatWidget(ToolOptionFloat* opts)
{
  QLineEdit* line=new QLineEdit(this);
  QDoubleValidator* vali=new QDoubleValidator(line);
  vali->setRange(opts->GetLowerLimit(), opts->GetUpperLimit(), 2); //allow two decimal places
  line->setValidator(vali);
  line->setText(QString::number(opts->GetValue()));
  connect(line, SIGNAL(textChanged(const QString&)), 
          new FloatOptBinder(line, opts), SLOT(TextChanged(const QString&)));  
  return line;
}

QWidget* ToolOptionsWidget::MakeButtonWidget(ToolOptionButton* opts)
{
  QPushButton* button=new QPushButton();
  button->setText(QString(opts->GetVerboseName().c_str()));
  connect(button, SIGNAL(released()), opts->GetReceiver(), opts->GetSlotMethod());
  return button;
}

void ToolOptionsWidget::Populate()
{
  QGridLayout* grid=new QGridLayout(this);
  grid->setSpacing(5);
  const ToolOptionGroups& groups=options_->GetGroups();
  int current_row=0;
  for (ToolOptionGroups::const_iterator i=groups.begin(), 
       e=groups.end(); i!=e; ++i) {
    ToolOptionGroup* group=*i;
    if (group->GetVerboseName()!="") {
      QLabel* group_header=new QLabel(QString(group->GetVerboseName().c_str()),
                                      this);
      group_header->setBackgroundRole(QPalette::AlternateBase);
      grid->addWidget(group_header, current_row++, 0, 1, 2);
    }
    for (ToolOptionList::const_iterator j=group->GetOptions().begin(),
         e2=group->GetOptions().end(); j!=e2; ++j) {
      ToolOption* opt=*j;
      QLabel* opt_label=new QLabel(QString(opt->GetVerboseName().c_str()),
                                   this);
      QWidget* widget=NULL;
      switch(opt->GetType()) {
        case ToolOption::INT:
          widget=MakeIntWidget(dynamic_cast<ToolOptionInt*>(opt));
          break;
        case ToolOption::FLOAT:
          widget=MakeFloatWidget(dynamic_cast<ToolOptionFloat*>(opt));
          break;
        case ToolOption::ENUM:
          widget=MakeEnumWidget(dynamic_cast<ToolOptionEnum*>(opt));
          break;
        case ToolOption::BUTTON:
          widget=MakeButtonWidget(dynamic_cast<ToolOptionButton*>(opt));
          break;
      }
#if defined(__APPLE__)
      opt_label->setAttribute(Qt::WA_MacSmallSize, true);
      widget->setAttribute(Qt::WA_MacSmallSize, true);  
#endif
      grid->addWidget(opt_label, current_row, 0, 1, 1,
                      Qt::AlignLeft|Qt::AlignTop);
      grid->addWidget(widget, current_row, 1, 1, 1, Qt::AlignLeft|Qt::AlignTop);
      
      current_row++;
      opt_label->setBuddy(widget);
    }
  }
  grid->addItem(new QSpacerItem(1, 1, QSizePolicy::Minimum, 
                                QSizePolicy::Expanding), current_row, 0, 1, 2);
}

}}
