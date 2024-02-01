//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  Authors: Ansgar Philippsen, Andreas Schenk
*/
#include "data_viewer_panel.hh"
#include "pointlist_overlay_base.hh"
#include "strategies.hh"

#include <QActionGroup>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QDialogButtonBox>
#include <QGroupBox>
#include <QSpinBox>
#include <QPixmap>
#include <QIcon>
#include <QColorDialog>


namespace ost { namespace img { namespace gui {

PointlistOverlayBaseSettings::PointlistOverlayBaseSettings(const QColor& ac, const QColor& pc, int s1, int s2, QWidget* p):
  QDialog(p),
  active_color(ac),passive_color(pc),symbol_size(s1),symbol_strength(s2)
{
  QPixmap pm(32,32);
  pm.fill(ac);
  active_color_b_ = new QPushButton(QIcon(pm),"Active Color");
  connect(active_color_b_,SIGNAL(pressed()),this,SLOT(OnActiveColor()));

  pm.fill(pc);
  passive_color_b_ = new QPushButton(QIcon(pm),"Passive Color");
  connect(passive_color_b_,SIGNAL(pressed()),this,SLOT(OnPassiveColor()));

  QSpinBox* sb1 = new QSpinBox;
  sb1->setRange(1,99);
  sb1->setValue(s1);
  connect(sb1,SIGNAL(valueChanged(int)),this,SLOT(OnSymbolSize(int)));

  QSpinBox* sb2 = new QSpinBox;
  sb2->setRange(1,3);
  sb2->setValue(s2);
  connect(sb2,SIGNAL(valueChanged(int)),this,SLOT(OnSymbolStrength(int)));

  QHBoxLayout* hbl = new QHBoxLayout;
  hbl->addWidget(active_color_b_);
  hbl->addWidget(passive_color_b_);

  QGroupBox* gb = new QGroupBox("Symbol Colors",this);
  gb->setLayout(hbl);

  main_layout_ = new QVBoxLayout;
  main_layout_->addWidget(gb);

  hbl = new QHBoxLayout;
  hbl->addWidget(sb1);
  hbl->addWidget(sb2);

  gb = new QGroupBox("Symbol Size and Width",this);
  gb->setLayout(hbl);

  main_layout_->addWidget(gb);

  QDialogButtonBox* bbox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);

  connect(bbox, SIGNAL(accepted()), this, SLOT(accept()));
  connect(bbox, SIGNAL(rejected()), this, SLOT(reject()));

  QVBoxLayout* top = new QVBoxLayout;
  top->addLayout(main_layout_);
  top->addWidget(bbox);

  setLayout(top);
}

void PointlistOverlayBaseSettings::OnActiveColor()
{
  QColor col=QColorDialog::getColor(active_color); 
  if(col.isValid()) {
    active_color=col;
    QPixmap pm(32,32);
    pm.fill(col);
    active_color_b_->setIcon(QIcon(pm));
  }
}

void PointlistOverlayBaseSettings::OnPassiveColor()
{
  QColor col=QColorDialog::getColor(passive_color); 
  if(col.isValid()) {
    passive_color=col;
    QPixmap pm(32,32);
    pm.fill(col);
    passive_color_b_->setIcon(QIcon(pm));
  }
}

void PointlistOverlayBaseSettings::OnSymbolSize(int v)
{
  symbol_size=v;
}

void PointlistOverlayBaseSettings::OnSymbolStrength(int v)
{
  symbol_strength=v;
}

/////////////////////////////

PointlistOverlayBase::PointlistOverlayBase(const String& name):
  Overlay(name),
  active_color_(QColor(255,0,0)),
  passive_color_(QColor(100,100,100)),
  symbolsize_(5),
  symbolstrength_(1)
{
  strategy_=SymbolDrawingStrategyPtr(new SquareDrawingStrategy());

  QActionGroup* ag1 = new QActionGroup(this);
  ag1->setExclusive(true);
  a_sq_ = new QAction("S&quares",ag1);
  a_sq_->setCheckable(true);
  a_sq_->setChecked(true);
  a_ci_ = new QAction("&Circles",ag1);
  a_ci_->setCheckable(true);
  a_cr_ = new QAction("Display Crosshair",this);
  a_cr_->setCheckable(true);
  a_settings_ = new QAction("Settings...",this);

  menu_=new QMenu("Pointlist Overlay");
  menu_->addAction(a_sq_);
  menu_->addAction(a_ci_);
  menu_->addSeparator();
  menu_->addAction(a_cr_);
  menu_->addSeparator();
  menu_->addAction(a_settings_);
}

QMenu* PointlistOverlayBase::GetMenu()
{
  return menu_;
}

void PointlistOverlayBase::OnMenuEvent(QAction* e)
{
  if(e==a_sq_ || e==a_ci_ || e==a_cr_) {
    if(a_sq_->isChecked()) {
      if(a_cr_->isChecked()) {
        ReplaceStrategy<CrosshairSquareDrawingStrategy>();
      }else{
        ReplaceStrategy<SquareDrawingStrategy>();
      }
    } else {
      if(a_cr_->isChecked()) {
        ReplaceStrategy<CrosshairCircleDrawingStrategy>();
      }else{
        ReplaceStrategy<CircleDrawingStrategy>();
      }
    }
  } else if(e==a_settings_) {
    PointlistOverlayBaseSettings props(active_color_,passive_color_,symbolsize_,symbolstrength_,e->parentWidget());
    if(props.exec()==QDialog::Accepted) {
      SetProps(&props);
    }
  }
}

void PointlistOverlayBase::DrawPointList(QPainter& pnt, DataViewerPanel* dvp, 
                                         const QColor& color,
                                         const std::vector<QPoint>& pointlist)
{
  strategy_->SetPenColor(color);
  strategy_->SetSymbolStrength(symbolstrength_);
  strategy_->SetSymbolSize(static_cast<int>((GetSymbolSize()-0.5) * dvp->GetZoomScale()));
  for(unsigned int i = 0; i < pointlist.size(); i++) {
    strategy_->Draw(pnt,pointlist[i]);
  }
  strategy_->SetSymbolSize(static_cast<int>((GetSymbolSize()-0.5) * dvp->GetZoomScale()));
}
void PointlistOverlayBase::DrawVariableSizePointList(QPainter& pnt, DataViewerPanel* dvp, 
                                                     const QColor& color,
                                                     const std::vector<std::pair<QPoint,double> >& pointlist)
{
  strategy_->SetPenColor(color);
  strategy_->SetSymbolStrength(symbolstrength_);
  for(unsigned int i = 0; i < pointlist.size(); i++) {
    strategy_->SetSymbolSize(static_cast<int>((GetSymbolSize()-0.5) * dvp->GetZoomScale()*pointlist[i].second));
    strategy_->Draw(pnt,pointlist[i].first);
  }
  strategy_->SetSymbolSize(static_cast<int>((GetSymbolSize()-0.5) * dvp->GetZoomScale()));
}

bool PointlistOverlayBase::GetCrosshair() const
{
  return a_cr_->isChecked();
}

void PointlistOverlayBase::SetCrosshair(bool flag)
{
  if(a_cr_->isChecked())	{
    if(a_sq_->isChecked()){
      ReplaceStrategy<CrosshairSquareDrawingStrategy>();
    }else{
      ReplaceStrategy<CrosshairCircleDrawingStrategy>();
    }
  } else {
    if(a_sq_->isChecked()){
      ReplaceStrategy<SquareDrawingStrategy>();
    }else{
      ReplaceStrategy<CircleDrawingStrategy>();
    }
  }
}

unsigned int PointlistOverlayBase::GetSymbolSize() const 
{
  return symbolsize_;
}

void PointlistOverlayBase::SetSymbolSize(unsigned int symbolsize)
{
  symbolsize_=symbolsize;
}

unsigned int PointlistOverlayBase::GetSymbolShape() const
{
  return a_sq_->isChecked() ? 1 : 0;
}

void PointlistOverlayBase::SetSymbolShape(unsigned int symbolshape)
{
  switch(symbolshape){
  case 1:
    if(a_cr_->isChecked()){
      ReplaceStrategy<CrosshairSquareDrawingStrategy>();
    }else{
      ReplaceStrategy<SquareDrawingStrategy>();
    }
    break;
  case 0:
    if(a_cr_->isChecked()){
        ReplaceStrategy<CrosshairCircleDrawingStrategy>();
    }else{
      ReplaceStrategy<CircleDrawingStrategy>();
    }
    break;
  }
  
}

template <class StrategyClass>
void PointlistOverlayBase::ReplaceStrategy()
{
  SymbolDrawingStrategyPtr s(new StrategyClass(*strategy_));
  strategy_.swap(s);
}

void PointlistOverlayBase::SetProps(PointlistOverlayBaseSettings* props)
{
  active_color_=props->active_color;
  passive_color_=props->passive_color;
  symbolsize_=props->symbol_size;
  symbolstrength_=props->symbol_strength;
}

void PointlistOverlayBase::SetActiveColor(const QColor& col)
{
  active_color_=col;
}

void PointlistOverlayBase::SetPassiveColor(const QColor& col)
{
  passive_color_=col;
}

}}}  //ns
