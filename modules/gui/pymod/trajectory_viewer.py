"""
This module defines a widget used for visualization of trajectories.

"""

__all__=('TrajWidget','SetTimeUnit','AddTrajectory','RemoveTrajectory','SetSpeed','ReleaseGfxEntity','FixGfxEntity')

#-------------------------------------------------------------------------------
# This file is part of the OpenStructure project <www.openstructure.org>
#
# Copyright (C) 2008-2020 by the OpenStructure authors
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#-------------------------------------------------------------------------------

import PyQt5 as _PyQt5
from PyQt5.QtCore import Qt as _Qt
from PyQt5.QtWidgets import QWidget as _QWidget
from PyQt5.QtWidgets import QLabel as _QLabel
from PyQt5.QtWidgets import QPushButton as _QPushButton
import ost as _ost
import math as _math


class TrajWidget(_QWidget):
  def __init__(self, trajlist=None, golist=None, parent=None,ref_index=0):
    _QWidget.__init__(self, parent, _Qt.Tool)
    self.setFocusPolicy(_Qt.ClickFocus)
    self.ref_index_=ref_index
    self.trajlist_=trajlist
    self.golist_=golist
    self.ehlist_=[go.GetView().handle for go in self.golist_]
    self.index_dict={}
    self.fix_dict={}
    self.modifiers=None
    for i,go in enumerate(self.golist_):
      self.index_dict[go.name]=i
    vb=_PyQt5.QtWidgets.QVBoxLayout()
    hb=_PyQt5.QtWidgets.QHBoxLayout()
    hb1=_PyQt5.QtWidgets.QHBoxLayout()
    hb2=_PyQt5.QtWidgets.QHBoxLayout()
    hb3=_PyQt5.QtWidgets.QHBoxLayout()
    hb4=_PyQt5.QtWidgets.QHBoxLayout()
    self.callback=None   
    self._slider=_PyQt5.QtWidgets.QSlider(self)
    self._slider.setOrientation(_Qt.Horizontal)
    self._speed_slider=_PyQt5.QtWidgets.QSlider(self)
    self._speed_slider.setOrientation(_Qt.Horizontal)
    self._speedLabel=_QLabel(self)
    self._speedLabel.setText('Speed:')
    self._speedLabel.setAlignment(_Qt.AlignLeft)
    self._play=_PyQt5.QtWidgets.QToolButton(self)
    self._repeat=_PyQt5.QtWidgets.QCheckBox(self)
    self._frame=_QLabel(self)
    self._frameNo=_QLabel(self)
    self._frameEnd=_QLabel(self)
    self._repeat.setText('Repeat')
    self._slider.setTracking(True)
    self._play.setText('Play')
    self._play.setCheckable(True)
    self._frame.setText('Frame: ')
    self._frameNo.setNum(0)
    self._frameNo.setAlignment(_Qt.AlignRight)
    self._frameEnd.setAlignment(_Qt.AlignLeft)
    self._timeLabel=_QLabel(self)
    self._timeNo=_QLabel(self)
    self._timeEnd=_QLabel(self)
    self._timeUnit=_QLabel(self)
    self._timeLabel.setText('Time: ')
    self._timeNo.setNum(0)
    self._timeNo.setAlignment(_Qt.AlignRight)
    self._timeEnd.setAlignment(_Qt.AlignLeft)
    self._timeUnit.setAlignment(_Qt.AlignLeft)
    self._speed_slider.setTracking(True)
    self._speed_slider_min=-50
    self._speed_slider_max=-10
    self._speed_slider.setRange(self._speed_slider_min,self._speed_slider_max)
    self._right_arrow=_QPushButton(">")
    self._right_arrow2=_QPushButton(">>")
    self._right_end=_QPushButton(">|")
    self._left_arrow=_QPushButton("<")
    self._left_arrow2=_QPushButton("<<")
    self._left_end=_QPushButton("|<")
    self._align_selection=_PyQt5.QtWidgets.QLineEdit()
    self._align=_QPushButton("Align")
    self._align_label=_QLabel(self)
    self._align_label.setText('Selection:')
    self._ref_entity_selection=_PyQt5.QtWidgets.QComboBox()
    for go in self.golist_:self._ref_entity_selection.addItem(go.name)
    self._ref_entity_selection_label=_QLabel(self)
    self._ref_entity_selection_label.setText('Ref. Entity:')
    hb.addWidget(self._play)
    hb.addWidget(self._repeat)
    hb.addWidget(self._ref_entity_selection_label)
    hb.addWidget(self._ref_entity_selection)
    hb1.addWidget(self._left_end)
    hb1.addWidget(self._left_arrow2)
    hb1.addWidget(self._left_arrow)
    hb1.addWidget(self._right_arrow)
    hb1.addWidget(self._right_arrow2)
    hb1.addWidget(self._right_end)
    hb2.addWidget(self._frame)
    hb2.addWidget(self._frameNo)
    hb2.addWidget(self._frameEnd)
    hb2.addWidget(self._timeLabel)
    hb2.addWidget(self._timeNo)
    hb2.addWidget(self._timeEnd)
    hb2.addWidget(self._timeUnit)
    hb3.addWidget(self._speedLabel)
    hb3.addWidget(self._speed_slider)
    hb4.addWidget(self._align_label)
    hb4.addWidget(self._align_selection)
    hb4.addWidget(self._align)
    self.setLayout(vb)
    vb.addLayout(hb)
    vb.addLayout(hb1)
    vb.addWidget(self._slider)
    vb.addLayout(hb2)
    vb.addLayout(hb3)
    vb.addLayout(hb4)
    self._SetSpeedSliderPos(self._speed_slider_min+0.5*(self._speed_slider_max-self._speed_slider_min))
    self._play.toggled.connect(self._TogglePlay)
    self._slider.valueChanged.connect(self._SliderValueChanged)
    self._speed_slider.valueChanged.connect(self._SpeedSliderValChanged)
    self._right_end.clicked.connect(self._RightEndClicked)
    self._right_arrow2.clicked.connect(self._Right2Clicked)
    self._right_arrow.clicked.connect(self._RightClicked)
    self._left_arrow2.clicked.connect(self._Left2Clicked)
    self._left_arrow.clicked.connect(self._LeftClicked)
    self._left_end.clicked.connect(self._LeftEndClicked)
    self._align.clicked.connect(self._AlignClicked)
    self._ref_entity_selection.currentIndexChanged.connect(self._EntitySelected)
    self._slider.setMinimum(0)
    self.frame_number_=self.traj_.GetFrameCount()-1
    self.timestep_=self.traj_.GetDelta()
    self.SetTimeUnit('ns')
    self.SetReferenceIndex(0)
    
  def _p2u(self,u):
    if u=="s":
      return 1.0e-12
    elif u=="ms":
      return 1.0e-9
    elif u=="us":
      return 1.0e-6
    elif u=="ns":
      return 1.0e-3
    elif u=="ps":
      return 1.0
    elif u=="fs":
      return 1.0e3
    raise RuntimeError("expected one of s,ms,us,ns,ps or fs for unit")

  def SetTimeUnit(self,u):
    self._time_prefactor=self._p2u(u)
    self._timeEnd.setText('/ '+ str((self.frame_number_-1)*self.timestep_*self._time_prefactor))
    self._timeNo.setNum(self.current_frame*self.timestep_*self._time_prefactor)
    self._timeUnit.setText('['+u+']')
  
  def _SetSpeedSliderPos(self,pos):
    self._speed_slider.setSliderPosition(pos)
    self._SpeedSliderValChanged(pos)
    
  def _SpeedSliderValChanged(self,speed_pos):
    self.time=_math.exp(-0.15*speed_pos)
    if self._play.isChecked():
      self._TogglePlay(False)
      self._TogglePlay(True)
  
  def _SetTime(self,t):
    self.time=t
    self._speed_slider.setSliderPosition(-1./0.15*_math.log(t))
    if self._play.isChecked():
      self._TogglePlay(False)
      self._TogglePlay(True)
  
  
  def _SliderValueChanged(self, pos):
    self.current_frame=pos
    for traj,go in zip(self.trajlist_,self.golist_):
      if go.name in self.fix_dict:continue
      traj.CopyFrame(self.current_frame)
      go.UpdatePositions()

  def _GetCurrentFrame(self):
    return self._slider.sliderPosition()
    
  def _SetCurrentFrame(self, pos):
    if self._slider.maximum()<pos:
      if self._repeat.isChecked():
        pos=0
      else:
        pos=self._slider.maximum()
    self._slider.setSliderPosition(pos)
    self._frameNo.setNum(pos)
    self._timeNo.setNum(pos*self.timestep_*self._time_prefactor)
    
  current_frame=property(_GetCurrentFrame, _SetCurrentFrame)
  
  def _GetReferenceTraj(self):
    return self.trajlist_[self.ref_index_]
  traj_=property(_GetReferenceTraj)
  
  def timerEvent(self, event):
    #if self.callback:
    #  self.callback(self.golist_)
    self.current_frame+=1
    for go in self.golist_:
      if go.name in self.fix_dict:continue
      go.BlurSnapshot()
      go.UpdatePositions()
  
  def _TogglePlay(self, playing):
    if playing:
      self.timer_id_=self.startTimer(self.time)
    else:
      self.killTimer(self.timer_id_)
  
  def _LeftClicked(self):
    if self.current_frame>0:
      self.current_frame-=1
  
  def _RightClicked(self):
    if self.current_frame<self.frame_number_-1:
      self.current_frame+=1
  
  def _Left2Clicked(self):
    if self.current_frame>=10:
      self.current_frame-=10
  
  def _Right2Clicked(self):
    if self.current_frame<self.frame_number_-10:
      self.current_frame+=10
  
  def _LeftEndClicked(self):
    if self.current_frame>0:
      self.current_frame=0
  
  def _RightEndClicked(self):
    if self.current_frame<self.frame_number_-1:
      self.current_frame=self.frame_number_-1
    
  def _AlignClicked(self):
    ref_eh=self.ehlist_[self.ref_index_]
    ref_v=ref_eh.Select(str(self._align_selection.text()))
    if ref_v.GetAtomCount()<=3:
      print('not enough atoms for alignment')
      return
    for i,t,eh in zip(list(range(len(self.trajlist_))),self.trajlist_,self.ehlist_):
      t=_ost.mol.alg.SuperposeFrames(t,eh.Select(str(self._align_selection.text())),ref_v)
      self.trajlist_[i]=t
      
  def _SetBlur(self, blur):
    for go in self.golist_:
      go.SetBlur(blur)

  def _GetBlur(self):
    return self.gfx_entity.GetBlur()

  blur=property(_GetBlur, _SetBlur)

  def SetReferenceIndex(self,ref_index):
    if type(ref_index)==type(''):
      self.ref_index_=self.index_dict[ref_index]
    else:
      self.ref_index_=ref_index
    self.frame_number_=self.traj_.GetFrameCount()
    self.timestep_=self.traj_.GetDelta()
    self._slider.setMaximum(self.frame_number_-1)
    self._frameEnd.setText('/ '+ str(self.frame_number_-1))
    self._timeEnd.setText('/ '+ str((self.frame_number_-1)*self.timestep_*self._time_prefactor))
    self._ref_entity_selection.setCurrentIndex(self.ref_index_)
  
  def _EntitySelected(self,index):
    ref_index=self.index_dict[str(self._ref_entity_selection.currentText())]
    self.SetReferenceIndex(ref_index)
  
  def keyPressEvent(self, event):
    key=event.key()
    if event.modifiers()==_Qt.ControlModifier:
      self.modifiers=event.modifiers()
    if self.modifiers==_Qt.ControlModifier:
      if key==_Qt.Key_Left:
        self._Left2Clicked()
      elif key==_Qt.Key_Right:
        self._Right2Clicked()
    else:
      if key==_Qt.Key_Left:
        self._LeftClicked()
      elif key==_Qt.Key_Right:
        self._RightClicked()
      elif key==_Qt.Key_Space:
        self._play.setChecked(not self._play.isChecked())
      else:
        _QWidget.keyPressEvent(self, event)
  
  def keyReleaseEvent(self, event):
    if event.key()==_Qt.Key_Control:
      self.modifiers=None
    else:
      _QWidget.keyPressEvent(self, event)
  
  def FixGfxEntity(self,index,frame_number=None):
    if not frame_number:frame_number=self.current_frame
    if type(index)==type(''):
      self.fix_dict[index]=frame_number
      i=self.index_dict[index]
    else:
      i=index
      self.fix_dict[self.golist_[i].name]=frame_number
    self.trajlist_[i].CopyFrame(frame_number)
    self.golist_[i].UpdatePositions()
  
  def ReleaseGfxEntity(self,index):
    if type(index)==type(''):self.fix_dict.pop(index)
    else:self.fix_dict.pop(self.golist_[index].name)
  
  def AddTrajectory(self,traj,go):
    self.trajlist_.append(traj)
    self.golist_.append(go)
    self.ehlist_.append(go.GetView().handle)
    self.index_dict[go.name]=len(self.golist_)-1
  
  def RemoveTrajectory(self,index):
    if type(index)==type(''):
      index=self.index_dict[index]
    self.trajlist_.pop(index)
    self.golist_.pop(index)
    self.ehlist_.pop(index)
  
  def SetSpeed(self,val):
  #Value should be between 0 and 1
    if not (val<=1. and val >=0.):
      print('Speed should be set between 0 and 1')
      return
    else:
      val=self._speed_slider_min-val*(self._speed_slider_min-self._speed_slider_max)
      self._SetSpeedSliderPos(val)
      print(val)
      
