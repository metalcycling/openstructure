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
  Authors: Ansgar Philippsen, Johan Hebert, Andreas Schenk
*/

#ifndef DATA_VIEWER_PANEL_BASE_HH_
#define DATA_VIEWER_PANEL_BASE_HH_

#include <map>    

#ifndef Q_MOC_RUN

#include <boost/shared_ptr.hpp>
#include <ost/base.hh>
#include <ost/img/extent.hh>
#include <ost/img/data_observer.hh>
#include <ost/img/normalizer_impl.hh>
#include <ost/img/raster_image/raster_image.hh>

#include <ost/gui/module_config.hh>
#include <ost/gui/data_viewer/viewer_normalizer.hh>

#include <ost/img/normalizer_impl.hh>

#endif

#include <QWidget>
#include <QCursor>
#include <QMenu>
#include <QPen>
#include <QBrush>
#include <QRubberBand>

namespace ost { namespace img { namespace gui {

    
class DLLEXPORT_OST_GUI DataViewerPanelBase: public QWidget, 
                                             public DataObserver 
{
  Q_OBJECT;

public:
  DataViewerPanelBase(const ImageHandle& data,QWidget* parent);
  virtual ~DataViewerPanelBase();

  virtual void SetData(const ImageHandle& d);

  //////////////////
  // event handling
  virtual void paintEvent(QPaintEvent* event);
  virtual void resizeEvent(QResizeEvent* event);
  virtual void wheelEvent (QWheelEvent* event);
  virtual void mousePressEvent(QMouseEvent* event);
  virtual void mouseReleaseEvent(QMouseEvent* event);
  virtual void mouseMoveEvent(QMouseEvent* event);
  virtual void keyPressEvent(QKeyEvent * event);


  //////////////////////
  // observer interface
  virtual void ObserverUpdate();
  virtual void ObserverUpdate(const Extent& e);
  virtual void ObserverUpdate(const Point& p);
  virtual void ObserverRelease();

  /////////////////
  // other methods

  //! update view
  void UpdateView(bool update_raster_image=true);
  void UpdateView(const QRect& rect,bool update_raster_image=true);

  //! retrieve ptr to internal normalizer
  ViewerNormalizerPtr GetNormalizer() const;

  //! do a renormalization with the current normalizer
  void Renormalize();
  

  //! return currently active selection
  /*!
    if no selection is present, return Extent()
  */
  Extent GetSelection() const;

  //! set currently active selection
  void SetSelection(const Extent& extent);

  //! convert window coordinates to image point
  Point WinToPoint(int mx, int my) const;
  Point WinToPoint(const QPoint& p) const;

  //! convert window coordinates to fractional image point
  geom::Vec2 WinToFracPoint(const QPoint& p) const;
  geom::Vec2 WinToFracPoint(int mx, int my) const;
  geom::Vec2 WinToFracPointCenter(const QPoint& p) const;
  geom::Vec2 WinToFracPointCenter(int mx, int my) const;

  //! convert image point to window coordinates
  QPoint PointToWin(const Point& p) const;
  //! convert fractional image point to window coordinates (top left corner of point)
  QPoint FracPointToWin(const geom::Vec2& p) const;
  //! convert fractional image point to window coordinates (center of point)
  QPoint FracPointToWinCenter(const geom::Vec2& p) const;

  //! checks if given point is visible in current panel
  bool IsWithin(const QPoint& p) const;

  //! retrieve the zoomscale
  Real GetZoomScale() const;
    
  //! set the zoomscale
  void SetZoomScale(Real zoom);
    
  //! retrieve the extent of the image
  Extent GetExtent();

  //! retrieve the pixel sampling of the image
  geom::Vec3 GetPixelSampling(); 

  
  //! set default cursor
  void SetDefaultCursor(const QCursor& c);
  void SetDefaultCursor(int id);

  //! set cursor
  void SetCursor(const QCursor& c);
  void SetCursor(int id);

  //! set display pixelvalues
  void SetDisplayPixelValues(bool show);
  bool GetDisplayPixelValues(){return display_pixel_values_;};
 
  //! re-center with spatial origin in the middle of the window
  void Recenter();

  
  /////////////////
  // other methods
  //! set coloring mode to either greyscale or phasecolor
  void SetColorMode(RasterImage::Mode m);
  RasterImage::Mode GetColorMode(){return cmode_;};

  //! set/get fast drawing modes for low and high magnification
  void SetFastLowMagnificationDrawing(bool flag);
  bool GetFastLowMagnificationDrawing();
  void SetFastHighMagnificationDrawing(bool flag);
  bool GetFastHighMagnificationDrawing();

  geom::Vec2 GetMousePosition(){return mouseposition_;};
  void MoveTo(const geom::Vec2& p);
  ImageHandle Extract(const Extent& e);
  void HideClickedPosition();
  void ShowClickedPosition(bool show);
  bool ClickedPositionIsVisible();
  bool HasSelection();
  geom::Vec3 GetClickedPosition();
  Real GetDataMin() const;
  Real GetDataMax() const;
  bool GetInvert() const;
  void SetInvert(bool invert);
  Real GetGamma() const;
  void SetGamma(Real gamma);
  void SetViewerMin(Real min);
  Real GetViewerMin() const;
  void SetViewerMax(Real max);
  Real GetViewerMax() const;
  void UpdateNormalizer(Real min, Real max, Real gamma, bool invert);
  int GetSlab();
  void SetSlab(int slab);
  void SetSelectionMode(int mode);
  int GetSelectionMode();
  void SetAntialiasing(bool f);
  bool GetAntialiasing() const;
  geom::Vec2 GetOffset() const;
  void SetOffset(const geom::Vec2& offset);



signals:
  void clicked(const geom::Vec3& mousepos);
  void zoomed(int zoomlevel);
  void selected(const Extent& selection);
  void deselected();
  void slabChanged(int slab);
  void released(void);

protected:
  QMenu* popupmenu_;

protected slots:
  // slots for popup menu
  void toggle_phase_color();
  void normalize();
  void center();
  void toggle_invert();
  void toggle_display_pixels();
  void toggle_show_click_pos();
  void zoom_plus();
  void zoom_minus();
  void zoom_reset();
  void slab_plus();
  void slab_minus();
  void slab_start();
  void slab_end();
  void slab_center();
  void toggle_fast_low_mag();
  void toggle_fast_high_mag();
  void toggle_antialiasing();


private:
  Real data_min_;
  Real data_max_;
  Point data_min_pos_;
  Point data_max_pos_;
  ViewerNormalizerPtr normalizer_;
  QImage* image_;
  QPixmap* pixmap_;
  QRubberBand* rubberband_;
  QPoint lastmouse_;
  int zoom_level_;
  bool update_raster_image_;
  Real offset_x_, offset_y_;

  geom::Vec3 clicked_position_; 
  bool show_clicked_position_;

  Real zoom_scale_;
  Real i_zoom_scale_;
  int slab_;

  QCursor cursor_;
  bool display_pixel_values_;
  geom::Vec2 mouseposition_;

  int last_x_,last_y_;
  int right_press_x_,right_press_y_;
  Extent selection_;
  int selection_mode_;
  
  RasterImage::Mode cmode_;
  std::map<int, QCursor> cursors_;

  bool use_update_extent_;
  Extent update_extent_;
  bool fast_low_mag_;
  bool fast_high_mag_;

  bool antialiasing_;
  Point drag_start_;
  QAction* phase_toggle_action_;
  QAction* invert_toggle_action_;
  QAction* show_click_pos_toggle_action_;
  QAction* fast_low_mag_toggle_action_;
  QAction* fast_high_mag_toggle_action_;
  QAction* display_pixel_toggle_action_;
  QAction* antialiasing_toggle_action_;


  void move(int dx, int dy);
  void slab(int dz);
  void zoom(int d);
  void extract_ri();

  void draw_extent(QPainter& p);
  void draw_pixel_values(QPainter& p);
  void update_min_max();
  void update_rubberband_from_selection_();
};

}}}  //ns

#endif

