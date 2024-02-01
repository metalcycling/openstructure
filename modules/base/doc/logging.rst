Logging
================================================================================

.. currentmodule:: ost

OpenStructure has a logging system going beyond what print statements can offer. Messages can be logged to the terminal to the graphical user interface or a file. Depending on the needs, groups of messages can be turned and off.

.. function:: LogError(message)
              LogWarning(message)
              LogScript(message)
              LogInfo(message)
              LogVerbose(message)
              LogDebug(message)
              LogTrace(message)
  
  Add a message to the log. For the choice of the appropriate logging level, 
  see :ref:`picking-logging-level`.
  
  :param message: The message to be logged
  :type  message: convertible to string
  
  Similar to the built-int print function, several arguments can be passed to 
  the logging functions. The arguments will be converted to string and then 
  concatenated together, separated by space. For example:
  
  .. code-block:: python
    
    LogMessage('here be numbers', 1, 2)
    
  will be converted to 'here be numbers 1 2'.
  
.. note::

  In C++, the logging facility is implemented as a set of macros, called
  `LOG_ERROR`, `LOG_WARNING`, `LOG_SCRIPT`, `LOG_INFO`, `LOG_VERBOSE`,
  `LOG_DEBUG` and `LOG_TRACE`. The last two are only active when compiling with 
  debugging symbols. When debugging symbols are off, they expand to an empty 
  macro and thus don't create any overhead.


Verbosity Level
--------------------------------------------------------------------------------

Several verbosity levels are available. Verbosity levels are represented
by an enumeration of integer values. They are wrapped to objects with
memorable names by the :class:`LogLevel` class. The available levels are
are summarized in the table below.

===========   ================   ===================
Level name    Verbosity value    LogLevel object
===========   ================   ===================
Error         0                  LogLevel.Error
Warning       1                  LogLevel.Warning
Script        2                  LogLevel.Script
Info          3                  LogLevel.Info
Verbose       4                  LogLevel.Verbose
Debug         5                  LogLevel.Debug
Trace         6                  LogLevel.Trace
===========   ================   ===================

You can change the verbosity level  with the following two methods:

.. function:: PushVerbosityLevel(verbosity)

  Change the verbosity level to the given integer value or :class:`LogLevel` 
  enumeration object. All log events which
  have a severity above *verbosity* will not be shown. By default, the log level
  is 2, meaning that errors, warnings and script logging events are visible. By
  setting it to -1, you can disable all logging.

  :param verbosity: Desired verbosity level
  :type  verbosity: :class:`int`

  .. code-block:: python
    
       # Display warnings and errors:
       ost.PushVerbosityLevel(ost.LogLevel.Warning)

       # Disable all logging:
       ost.PushVerbosityLevel(-1)

.. function:: PopVerbosityLevel()

  Change the log level back to the previous verbosity level. It is an error to 
  pop the verbosity level without a matching call to 
  :func:`PushVerbosityLevel`.

.. function:: GetVerbosityLevel()

  :return: The current verbosity level
  :rtype:  :class:`int`

.. class:: LogLevel

  Enumerates the logging levels (see :ref:`picking-logging-level`). Values:

    .. attribute:: Error, Warning, Script, Info, Verbose, Debug, Trace

      The enumerated LogLevel object, which wraps the corresponding integer
      value. Note that these attributes are LogLevel objects themselves.

      :type: :class:`LogLevel`

      .. code-block:: python

          ost.LogLevel.Info
          # Outputs: ost._ost_base.LogLevel.Info

          print(ost.LogLevel.Info)
          # Outputs: Info
          
          int(ost.LogLevel.Info)
          # Outputs: 3

  These objects behave like integers, meaning that numeric comparisons work
  as expected. So for instance if you want to increase verbosity to the Info
  level, but leave it unchanged if it was already set to a higher value (such
  as Debug), you can do the following:

  .. code-block:: python

    new_level = max(ost.GetVerbosityLevel(), ost.LogLevel.Info)
    ost.PopVerbosityLevel()
    ost.PushVerbosityLevel(new_level)


         
Log sinks
--------------------------------------------------------------------------------

When running OpenStructure from the command-line, the log messages are by 
default output to stderr. When running DNG, the log messages are additionally 
logged to the messages widget. However, it is also possible to log into a file 
or theoretically even to a remote computer. All these are instances of so-called 
log sinks: classes that derive from LogSink and implement the LogMessage method.


.. class:: LogSink

  .. method:: LogMessage(message, severity)

     This method is called whenever something gets logged. This method must be 
     implemented by all subclasses.
     
     :param message: The logged message
    
     :type  message: :class:`str`
     
     :param severity: Marks how severe the logged message is. Errors have 
          severity 0, warnings 1 etc.
    
     :type  severity: :data:`int`

For convenience, there are 3 LogSink implementations available in OpenStructure that are sufficient for most use cases.

.. class:: FileLogSink(filename)

  The FileLogSink logs all messages into the given file.
  
  :param filename: The filename
  :type filename: :class:`str`

.. class:: StreamLogSink(stream)

  The stream log sink writes all log messages to the stream. stream must have a 
  write method that accepts a string. To write messages to stderr, use
  
  .. code-block:: python
  
    stderr_sink=ost.StreamLogSink(sys.stderr)
    ost.PushLogSink(stderr_sink)
    ost.LogInfo('Welcome, master')
  
.. class:: MultiLogSink

     A LogSink for multiplexing the log messages into multiple sinks at the same 
     time, e.g. the terminal and the messages widget.
     
  .. method:: AddSink(sink)
  
    Add a new sink. The sink's :meth:`LogSink.LogMessage` method will be called 
    every time something gets logged.
    
    :type sink: :class:`~ost.LogSink`
    :param sink: the log sink to be added

  .. method:: RemoveSink(sink)
  
    Remove the given sink. If the doesn't exist, this method has no effect.
    
    :type sink: :class:`~ost.LogSink`
    :param sink: the log sink to be removed
    

To change the current log sink you can use the following methods:

.. method:: PushLogSink(sink)

   Push the new sink onto the log sink stack. All of the messages will now be 
   logged to the new sink. To switch back to the previous log sink, use 
   :func:`PopLogSink`.

.. method:: PopLogSink()

   Change the log sink back to the previous one. It is an error to pop the log 
   sink when there is only one log sink on the stack.

.. method:: GetCurrentLogSink()

   Get the current (active) log sink.

.. _picking-logging-level:

Guidelines for picking logging level
--------------------------------------------------------------------------------

Each logging event has an associated level that marks its importance. For example, users should always see errors, but they do not need to see detailed information on the loading process. Here is a list of guidelines that we use in the code. We encourage developers to adhere to these guidelines as closely as possible.

**Error**:
  Very important message to the user, some command did not complete as expected
  or was aborted.

**Warning**:
  Diagnose potential problems that do not abort the execution, but may 
  point to a misconfiguration/misuse. This level is turned on by default.

**Script**:
  Logging level that should be used from scripts, e.g. to report progress. These 
  logging messages are turned on by default.

**Info**:
  Informative and important messages that summarize a complex command, such as 
  information on a loaded file, or results from an algorithm. These logging 
  messages are not turned on by default.

**Verbose**:
  Grey-zone between user and developer need, and perhaps the hardest to get
  right. This is the lowest logging level users will be able to see when they 
  use an optimized build. An example for this is the OpenGL setup/info in gfx, 
  or the path search during startup, or more detailed info on file IO. These 
  messages are not turned on by default.

**Debug**:
  For developers, but not quite at the trace level. This level is turned off by 
  default and only enabled when compiling with debugging symbols.

**Trace**:
  Used to debug inner loops. Once turned on, you will probably get more debug 
  output that you will be able to handle. This level is turned off by default
  and only enabled when compiling with debugging symbols.


Example
--------------------------------------------------------------------------------

The following snippet explains how to create a custom log sink which logs to the 
terminal (or the python shell in DNG). The logger also prints the current time.

.. code-block:: python

  import datetime
  class PyLogger(ost.LogSink):
    def __init__(self):
      ost.LogSink.__init__(self)

    def LogMessage(self, message, severity):
      levels=['ERROR', 'WARNING', 'SCRIPT', 'INFO', 
              'VERBOSE', 'DEBUG', 'TRACE']
      level=levels[severity]
      print('%s[%s]: %s' % (level, str(datetime.datetime.now()), message))

  py_logger=PyLogger()
  ost.PushLogSink(py_logger)

  ost.LogInfo("amazing logging system")
