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
 Author: Stefan Scheuber
 */

#include <ost/gui/gosty_app.hh>
#include <ost/platform.hh>

#include <ost/gui/widget_registry.hh>

#include "log_reader.hh"
#include "message_widget.hh"

#include <QApplication>
#include <QDir>
#include <QMenu>
#include <QVBoxLayout>

namespace ost {namespace gui {

class MessageWidgetFactory: public WidgetFactory {
public:
  MessageWidgetFactory() :
    WidgetFactory("ost::gui::MessageWidget", "Messages") {
  }

  virtual Widget* Create(QWidget* parent) {
    return GostyApp::Instance()->GetMessageWidget();
  }
};

OST_REGISTER_WIDGET(MessageWidget, MessageWidgetFactory);

MessageWidget::MessageWidget(QWidget* parent) :
  Widget(NULL, parent), model_(new QStandardItemModel(this)), view_(
      new QListView(this)) {
  QVBoxLayout* layout = new QVBoxLayout(this);
  layout->setMargin(0);
  layout->setSpacing(0);
  view_->setAttribute(Qt::WA_MacShowFocusRect, false);
  view_->setAttribute(Qt::WA_MacSmallSize, true);
  view_->setModel(model_);
  view_->setSelectionBehavior(QAbstractItemView::SelectRows);
  view_->setDragEnabled(true);
  view_->setContextMenuPolicy(Qt::CustomContextMenu);
  view_->setFrameShape(QFrame::NoFrame);
  view_->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);
  layout->addWidget(view_);

  connect(view_, SIGNAL(customContextMenuRequested(const QPoint&)), this,
      SLOT(ContextMenuRequested(const QPoint&)));

  QDir icon_path(GetSharedDataPath().c_str());
  icon_path.cd("gui");
  icon_path.cd("icons");

  QAction* clear_action = new QAction(this);
  clear_action->setToolTip("Clear info panel");
  clear_action->setIcon(QIcon(icon_path.absolutePath() + QDir::separator()
      + QString("delete_icon.png")));
  connect(clear_action, SIGNAL(triggered(bool)), this, SLOT(Clear()));
  this->actions_.append(clear_action);
  QObject::connect(&PythonInterpreter::Instance(),
                   SIGNAL(ErrorOutput(unsigned int, const QString &)),
                   this,
                   SLOT(ErrorOutput(unsigned int, const QString &)));
  new LogReader(this);
}

void MessageWidget::Update() {
  view_->viewport()->update();
}

void MessageWidget::LogMessage(const QString& message, QMessageBox::Icon icon) {
  QPixmap pix_icon = this->GetIcon(icon, this);
  QStandardItem* item = new QStandardItem();
  item->setData(QVariant(icon));
  item->setText(message);
  item->setIcon(QIcon(pix_icon));
  item->setEditable(false);
  this->model_->appendRow(item);
  this->Increase(icon);
}

void MessageWidget::ErrorOutput(unsigned int id,const QString& output)
{
  this->LogMessage(output, QMessageBox::Critical);
}

void MessageWidget::LogMessage(QStandardItem* item) {
  item->setData(QVariant(QMessageBox::NoIcon));
  this->model_->appendRow(item);
  this->Increase(QMessageBox::NoIcon);
}

void MessageWidget::LogMessage(const QString& message, QIcon icon) {
  QStandardItem* item = new QStandardItem();
  item->setText(message);
  item->setData(QVariant(QMessageBox::NoIcon));
  item->setIcon(icon);
  item->setEditable(false);
  this->model_->appendRow(item);
  this->Increase(QMessageBox::NoIcon);
}

void MessageWidget::Increase(QMessageBox::Icon icon){
  if(count_map_.contains(icon)){
    count_map_[icon] = count_map_[icon] + 1;
  }
  else{
    count_map_[icon] = 1;
  }
  emit MessageCountChanged(icon);
}

void MessageWidget::Decrease(QMessageBox::Icon icon){
  if(count_map_.contains(icon)){
    if(count_map_[icon]>0){
      count_map_[icon] = count_map_[icon] - 1;
    }
    else{
      count_map_[icon] = 0;
    }
  }
  emit MessageCountChanged(icon);
}

QPixmap MessageWidget::GetIcon(QMessageBox::Icon icon, QWidget* widget) {
  QDir icon_path(GetSharedDataPath().c_str());
  icon_path.cd("gui");
  icon_path.cd("icons");

  QStyle *style = widget ? widget->style() : QApplication::style();
  int icon_size = style->pixelMetric(QStyle::PM_MessageBoxIconSize, 0, widget);
  QIcon tmp_icon;
#if  defined(__APPLE__)
  switch (icon) {
    case QMessageBox::Information:
      tmp_icon = QIcon(icon_path.absolutePath()+QDir::separator()+QString("information_icon.png"));
      break;
    case QMessageBox::Warning:
      tmp_icon = QIcon(icon_path.absolutePath()+QDir::separator()+QString("warning_icon.png"));
      break;
    case QMessageBox::Critical:
      tmp_icon = QIcon(icon_path.absolutePath()+QDir::separator()+QString("critical_icon.png"));
      break;
    default:
      break;
  }
#else
  switch (icon) {
  case QMessageBox::Information:
    tmp_icon = style->standardIcon(QStyle::SP_MessageBoxInformation, 0, widget);
    break;
  case QMessageBox::Warning:
    tmp_icon = style->standardIcon(QStyle::SP_MessageBoxWarning, 0, widget);
    break;
  case QMessageBox::Critical:
    tmp_icon = style->standardIcon(QStyle::SP_MessageBoxCritical, 0, widget);
    break;
  case QMessageBox::Question:
    tmp_icon = style->standardIcon(QStyle::SP_MessageBoxQuestion, 0, widget);
  default:
    break;
  }
#endif
  if (!tmp_icon.isNull())
    return tmp_icon.pixmap(icon_size, icon_size);
  return QPixmap();
}

void MessageWidget::Clear() {
  this->model_->clear();
  this->count_map_.clear();
  emit AllCleared();
}

void MessageWidget::RemoveSelected() {
  QItemSelectionModel* selection_model = this->view_->selectionModel();
  const QItemSelection& item_selection = selection_model->selection();
  const QModelIndexList& model_indexes = item_selection.indexes();
  for (int i = 0; i < model_indexes.size(); i++) {
    int row = model_indexes[i].row();
    QMessageBox::Icon icon = QMessageBox::Icon(this->model_->item(row)->data().toInt());
    this->model_->removeRow(row);
    this->Decrease(icon);
  }
}

ActionList MessageWidget::GetActions() {
  return this->actions_;
}

void MessageWidget::ContextMenuRequested(const QPoint& pos) {
  QAction* remove_selected_action = new QAction("Remove", this);
  remove_selected_action->setToolTip("Remove this item");
  connect(remove_selected_action, SIGNAL(triggered(bool)), this,
      SLOT(RemoveSelected()));

  QMenu* menu = new QMenu();
  menu->addAction(remove_selected_action);
  if (! menu->actions().empty()) {
    menu->popup(view_->viewport()->mapToGlobal(pos));
  }
}

int MessageWidget::GetMessagesCount(QMessageBox::Icon icon) {
  if(count_map_.contains(icon)){
    return count_map_[icon];
  }
  return 0;
}

int MessageWidget::GetTotalMessagesCount() {
  return model_->rowCount();
}

MessageWidget::~MessageWidget() {}

}
} // ns
