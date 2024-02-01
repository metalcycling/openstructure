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
#include <boost/format.hpp>


#include <ost/img/image.hh>

#include "info_panel.hh"

#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QGroupBox>
namespace ost { namespace img { namespace gui {

InfoPanelLabel::InfoPanelLabel(const QString& t, QWidget* parent):
  QLabel(t,parent)
{
  setFrameStyle(QFrame::NoFrame);
  setContentsMargins(1,1,1,1);
  setIndent(0);
  setWordWrap(false);
  setMargin(0);
}

InfoPanel::InfoPanel(QWidget* parent):
  QWidget(parent),
  image_extent_(new InfoPanelLabel("",this)),
  image_size_(new InfoPanelLabel("",this)),
  image_sampling_(new InfoPanelLabel("",this)),
  image_type_(new InfoPanelLabel("",this)),
  mouse_xy_(new InfoPanelLabel("",this)),
  mouse_val_(new InfoPanelLabel("",this)),
  click_xy_(new InfoPanelLabel("",this)),
  click_val_(new InfoPanelLabel("",this)),
  selection_xy_(new InfoPanelLabel("",this)),
  size_xy_(new InfoPanelLabel("",this)),
  distance_(new InfoPanelLabel("",this)),
  popupmenu_(new QMenu(this)),
  amp_pha_(false),
  click_(geom::Vec2())
{
  QFont fnt=font();
  fnt.setPointSize(font().pointSize()-1);
  setFont(fnt);
  //int maxw = QFontMetrics(fnt).width("(9999,9999,1) to (9999,9999,1)");
  //setMaximumWidth(maxw);
  //setMinimumWidth(maxw);

  QFont fnt2(fnt);
  fnt2.setWeight(QFont::Bold);

  QVBoxLayout* top=new QVBoxLayout;
  top->setSpacing(0);
  QLabel* lab = new QLabel("Image");
  lab->setFont(fnt2);
  top->addWidget(lab);
  top->addWidget(image_extent_);
  top->addWidget(image_size_);
  top->addWidget(image_type_);
  top->addWidget(image_sampling_);

  lab = new QLabel("Mouse");
  lab->setFont(fnt2);
  top->addWidget(lab);
  top->addWidget(mouse_xy_);
  top->addWidget(mouse_val_);
  top->addWidget(click_xy_);
  top->addWidget(click_val_);
  top->addWidget(distance_);

  lab = new QLabel("Selection");
  top->addWidget(lab);
  lab->setFont(fnt2);
  top->addWidget(selection_xy_);
  top->addWidget(size_xy_);

  top->addStretch(1);

  setLayout(top);
}

void InfoPanel::UpdateMouseAndDistances(const geom::Vec2& point)
{
  int x=static_cast<int>(point[0]);
  int y=static_cast<int>(point[1]);
  int click_x=static_cast<int>(click_[0]);
  int click_y=static_cast<int>(click_[1]);
  std::ostringstream str;
  str << boost::format("%8.2f %8.2f") % point[0] % point[1];
  mouse_xy_->setText(QString::fromStdString(str.str()));
  if (click_ != geom::Vec2()) {
    Point point_distance= Point(x,y)-Point(click_x,click_y);
    geom::Vec2 position_distance= geom::Vec2(point[0]-click_[0],point[1]-click_[1]);
    Real dist1 = sqrt(static_cast<Real>(point_distance[0]*point_distance[0]+point_distance[1]*point_distance[1]));
    Real dist2 = sqrt(position_distance[0]*position_distance[0]+position_distance[1]*position_distance[1]);
    str.str("");
    str << boost::format("Distance:\n %8.2f %8.2f") % dist1 % dist2;
    distance_->setText(QString::fromStdString(str.str()));
  }
}

void InfoPanel::UpdateClick(const geom::Vec2& point)
{
  std::ostringstream str;
  str << boost::format("Clicked at: %8.2f %8.2f") % point[0] % point[1];
  click_xy_->setText(QString::fromStdString(str.str()));
  click_ = point;
  update();
}

void InfoPanel::SetSelection(const Extent& e)
{
  std::ostringstream str;
  Point sta=e.GetStart();
  Point end=e.GetEnd();
  str << boost::format("(%d,%d,%d) to (%d,%d,%d)") % sta[0] % sta[1] % sta[2] % end[0] % end[1] % end[2];
  selection_xy_->setText(QString::fromStdString(str.str()));
  str.str("");
  str << e.GetSize()[0] << "x" << e.GetSize()[1];
  if(e.GetSize()[2]>1) {
    str << "x" << e.GetSize()[2];
  }
  str << " pixel";
  size_xy_->setText(QString::fromStdString(str.str()));
  update();
}

void InfoPanel::SetMousePoint(const geom::Vec2& point, Real val)
{
  UpdateMouseAndDistances(point);
  std::ostringstream str;
  str << boost::format("Pixel Value: %8.2e") % val;
  mouse_val_->setText(QString::fromStdString(str.str()));
}

void InfoPanel::SetMousePoint(const geom::Vec2& point, Complex val)
{
  std::ostringstream str;
  UpdateMouseAndDistances(point);
  str << "Pixel Value: ";
  if(amp_pha_){
    str << boost::format("%8.2e , %8.2f") % std::abs(val) % std::arg(val);
  }else{
    Real im=imag(val);
    if(im>0){
      str << boost::format("%8.2e + %8.2ei") % real(val) % imag(val);
    }else{
      str << boost::format("%8.2e - %8.2ei") % real(val) % -imag(val);
    }
  }
  mouse_val_->setText(QString::fromStdString(str.str()));
  update();
}

void InfoPanel::SetClickPoint(const geom::Vec2& point, Complex val)
{
  std::ostringstream str;
  UpdateClick(point);
  str << "Pixel Value: ";
  if(amp_pha_){
    str << boost::format("%8.2e , %8.2f") % std::abs(val) % std::arg(val);
  }else{
    Real im=imag(val);
    if(im>0){
      str << boost::format("%8.2e + %8.2ei") % real(val) % imag(val);
    }else{
      str << boost::format("%8.2e - %8.2ei") % real(val) % -imag(val);
    }
  }
  click_val_->setText(QString::fromStdString(str.str()));
  update();
}

void InfoPanel::SetClickPoint(const geom::Vec2& point, Real val)
{
  UpdateClick(point);
  std::ostringstream str;
  str << boost::format("Pixel Value: %8.2e") % val;
  click_val_->setText(QString::fromStdString(str.str()));
  update();
}

void InfoPanel::ClearClickPosition()
{
  click_xy_->setText("");
  click_val_->setText("");
  distance_->setText("");
  update();
}

void InfoPanel::ClearSelection()
{
  selection_xy_->setText("");
  size_xy_->setText("");
  update();
}

void InfoPanel::SetImageInfo(const ConstData& data)
{
  std::ostringstream str;
  Extent ext=data.GetExtent();
  Point sta=ext.GetStart();
  Point end=ext.GetEnd();
  Size siz=ext.GetSize();
  str << boost::format("(%d,%d,%d) to (%d,%d,%d)") % sta[0] % sta[1] % sta[2] % end[0] % end[1] % end[2];
  image_extent_->setText(QString::fromStdString(str.str()));
  str.str("");
  str << siz[0] << "x" << siz[1];
  if(siz[2]>1) {
    str << "x" << siz[2];
  }
  str << " pixel";
  image_size_->setText(QString::fromStdString(str.str()));
  str.str("");
  str << boost::format("Sampling: %4.1eA") % data.GetPixelSampling()[0];
  image_sampling_->setText(QString::fromStdString(str.str()));
  str.str("");
  str << "Type: ";
  if(data.IsReal()) {str << "Real";} else {str<<"Complex";}
  str << " ";
  if(data.IsSpatial()) {str << "Spatial";} else {str<<"Frequency";}
  image_type_->setText(QString::fromStdString(str.str()));
}


}}}  //ns

