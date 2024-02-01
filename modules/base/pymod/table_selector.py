


class TableSelector:

  def __init__(self, col_types, col_names, query):

    self.col_types=col_types
    self.col_names=col_names
    self.query=query
    
    self.valid_operators=dict()

    self.valid_operators['and']={'type':'boolean','precedence':5}
    self.valid_operators['or']={'type':'boolean','precedence':6}

    self.valid_operators['!']={'type':'comparison','precedence':4,'extension':'!='}
    self.valid_operators['!=']={'type':'comparison','precedence':4,'extension':None}
    self.valid_operators['=']={'type':'comparison','precedence':4,'extension':'=='}
    self.valid_operators['==']={'type':'comparison','precedence':4,'extension':None}
    self.valid_operators['<=']={'type':'comparison','precedence':3,'extension':None}
    self.valid_operators['>=']={'type':'comparison','precedence':3,'extension':None}
    self.valid_operators['>']={'type':'comparison','precedence':3,'extension':'>='}
    self.valid_operators['<']={'type':'comparison','precedence':3,'extension':'<='}

    self.valid_operators['+']={'type':'arithmetic','precedence':2}
    self.valid_operators['-']={'type':'arithmetic','precedence':2}
    self.valid_operators['/']={'type':'arithmetic','precedence':1}
    self.valid_operators['*']={'type':'arithmetic','precedence':1}

    self.valid_operators['(']={'type':'left_bracket','precedence':float('NaN')}
    self.valid_operators['[']={'type':'left_bracket','precedence':float('NaN')}
    self.valid_operators['{']={'type':'left_bracket','precedence':float('NaN')}
    self.valid_operators[')']={'type':'right_bracket','precedence':float('NaN')}
    self.valid_operators[']']={'type':'right_bracket','precedence':float('NaN')}
    self.valid_operators['}']={'type':'right_bracket','precedence':float('NaN')}

    self.split_expression=self._ExpressionLexer(self.query)
    self.parsed_expression=self._ParseExpression(self.split_expression)
    self.rpn_expression=self._ShuntingYard(self.parsed_expression)

    self.tab_indices=list()
    self.exp_indices=list()

    #extract column indices
    for i, exp in enumerate(self.rpn_expression):
      if exp in self.col_names:
        self.tab_indices.append(self._GetIndex(exp))
        self.exp_indices.append(i)

  def EvaluateRow(self,row):
      for ti, ei in zip(self.tab_indices, self.exp_indices):
        #check for NaN
        if row[ti]!=row[ti]:
          self.rpn_expression[ei]=None
        else:
          self.rpn_expression[ei] = row[ti]
      if self._EvaluateRPN(list(self.rpn_expression)):
        return True
      return False

  def _GetIndex(self, col):
    if col not in self.col_names:
      raise ValueError('Table Selector has no column named "%s"' % col)
    return self.col_names.index(col)
    
  def _EvaluateAnd(self, lhs, rhs):
    return lhs==True and rhs==True

  def _EvaluateOr(self, lhs, rhs):
    return lhs==True or rhs==True

  def _EvaluateEqual(self, lhs, rhs):
    return lhs==rhs

  def _EvaluateNonEqual(self, lhs, rhs):
    return lhs!=rhs

  def _EvaluateLower(self, lhs, rhs):
    if lhs==None or rhs==None:
      return False
    return lhs<rhs

  def _EvaluateGreater(self, lhs, rhs):
    if lhs==None or rhs==None:
      return False
    return lhs>rhs

  def _EvaluateLowerEqual(self, lhs, rhs):
    if lhs==None or rhs==None:
      return False
    return lhs<=rhs

  def _EvaluateGreaterEqual(self, lhs, rhs):
    if lhs==None or rhs==None:
      return False
    return lhs>=rhs

  def _EvaluateAdd(self, lhs, rhs):
    if lhs==None or rhs==None:
      return None
    return lhs+rhs

  def _EvaluateSubtract(self, lhs, rhs):
    if lhs==None or rhs==None:
      return None
    return lhs-rhs

  def _EvaluateMultiply(self, lhs, rhs):
    if lhs==None or rhs==None:
      return None
    return lhs*rhs

  def _EvaluateDivide(self, lhs, rhs):
    if lhs==None or rhs==None:
      return None
    return lhs/rhs


  def _EvaluateOperator(self, op, lhs, rhs):

    #this function assumes, that all NaN values have been replaced by None!

    if op=='+':
      return self._EvaluateAdd(lhs, rhs)
    elif op=='-':
      return self._EvaluateSubtract(lhs, rhs)
    elif op=='/':
      return self._EvaluateDivide(lhs, rhs)
    elif op=='*':
      return self._EvaluateMultiply(lhs, rhs)
    elif op=='and':
      return self._EvaluateAnd(lhs, rhs)
    elif op=='or':
      return self._EvaluateOr(lhs, rhs)
    elif op=='=' or op=='==':
      return self._EvaluateEqual(lhs, rhs)
    elif op=='!=' or op=='!':
      return self._EvaluateNonEqual(lhs, rhs)
    elif op=='<':
      return self._EvaluateLower(lhs, rhs)
    elif op=='>':
      return self._EvaluateGreater(lhs, rhs)
    elif op=='<=':
      return self._EvaluateLowerEqual(lhs, rhs)
    elif op=='>=':
      return self._EvaluateGreaterEqual(lhs, rhs)

    else:
      raise ValueError('Unknown operator: '+op)

  def _EvaluateRPN(self, RPNExp):
    #Evaluates the reverse polish notation
    stack=list()
    while True:
      if len(RPNExp)==0:
        break
      exp=RPNExp.pop(0)
      if exp in self.valid_operators:
        if len(stack)<2:
          raise ValueError('Cannot evaluate operator on less than two operands!')
        rhs=stack.pop()
        lhs=stack.pop()
        result=self._EvaluateOperator(exp, lhs, rhs)
        if result==None:
          return False
        stack.append(result)
      else:
        stack.append(exp)
    if len(stack)>1:
      raise ValueError('Too many operands for given operators!')
    return stack.pop()

  def _ShuntingYard(self, split_expression):
    #Creates the so called reverse polish notation out of the expression parser output.
    #note, that there won't be parenthesis anymore and potential parenthesis
    #mismatches get recognized.
    #The shunting yard algorithm from dijkstra gets used.

    output_stack=list()
    operator_stack=list()

    while True:
      if len(split_expression)==0:
        while True:
          if len(operator_stack)==0:
            break
          if self.valid_operators[operator_stack[-1]]['type'] in ['left_bracket','right_bracket']:
            raise ValueError('Parenthesis mismatch!')
          output_stack.append(operator_stack.pop())
        break

      exp=split_expression.pop(0)

      if exp in self.valid_operators:
        if self.valid_operators[exp]['type']=='left_bracket':
          operator_stack.append(exp)
          continue

      if exp in self.valid_operators:
        if self.valid_operators[exp]['type'] == 'right_bracket':
          while True:
            if len(operator_stack)==0:
              raise ValueError('Parenthesis mismatch!')
            if self.valid_operators[operator_stack[-1]]['type']=='left_bracket':
              operator_stack.pop()
              break
            output_stack.append(operator_stack.pop())
          continue

      if exp in self.valid_operators:
        prec=self.valid_operators[exp]['precedence']
        while len(operator_stack)>0:
          if self.valid_operators[operator_stack[-1]]['type']=='left_bracket':
            break
          elif prec>=self.valid_operators[operator_stack[-1]]['precedence']:
            output_stack.append(operator_stack.pop())
          else:
            break
        operator_stack.append(exp)
        continue
      output_stack.append(exp)

    return output_stack

  def _ParseSubExpression(self, subexpression):

    valid_types={'float':'numeric','int':'numeric','string':'string','bool':'bool'}

    column_names=list()
    column_types=list()
   
    final_expression=list()


    for item in subexpression:
      if item in self.col_names:
        column_names.append(item)
        column_types.append(valid_types[self.col_types[self._GetIndex(item)]])
      
    unique_type=list(set(column_types))
    if len(unique_type)>1:
      raise ValueError('Try to compare columns '+','.join(column_names)+' which have inconsistent types!')
    if len(unique_type)==0:
      raise ValueError('Try to evaluate subexpression '+' '.join(subexpression)+' that contains no valid column name of current table!')
    
    for item in subexpression:
      if item in self.valid_operators:
        final_expression.append(item)
        continue
      if item in column_names:
        final_expression.append(item)
        continue
      if unique_type[0]=='numeric':
        if item in ['NaN','nan','None','none']:
          final_expression.append(None)
          continue
        else:
          try:
            final_expression.append(float(item))
            continue
          except:
            raise RuntimeError('Tried to cast '+item+' into numeric type to compare with column(s) '+','.join(column_names)+', but failed!')
      elif unique_type[0]=='bool':
        if item in ['None','none']:
          final_expression.append(None)
          continue
        if item in ['true','True']:
          final_expression.append(True)
          continue
        if item in ['false','False']:
          final_expression.append(False)
          continue
        raise RuntimeError('Tried to cast '+item+' into boolean type to compare with column(s) '+','.join(column_names)+', but failed!')
      elif unique_type[0]=='string':
        final_expression.append(item)

    return final_expression


  def _ParseExpression(self, split_expression):

    #check for problematic cases like 'a<=b<=c'. We don't know which operator to evaluate first
    for i in range(len(split_expression)-3):
      if (split_expression[i] in self.valid_operators) and (split_expression[i+2] in self.valid_operators):
        if self.valid_operators[split_expression[i]]['precedence']==self.valid_operators[split_expression[i+2]]['precedence']:
          raise ValueError('Cannot Evaluate '+' '.join(split_expression[i:i+3])+' since both operators have same precedence!')

    #handle , operator
    #replaces an expression like 'rnum=1,2,3' with '(rnum=1 or rnum=2 or rnum=3)'

    temp_split_expression=list()
    skips=0

    for i in range(len(split_expression)):
      if skips>0:
        skips-=1
        continue
      if ',' in split_expression[i]:

        if split_expression[max(0,i-1)] != '=' and split_expression[min(i+1,len(split_expression)-1)] != '=':
          raise ValueError('Can evaluate \',\' operator only in combination with \"=\" in subexpression ',' '.join(split_expression[max(0,i-1):min(i+1,len(split_expression))]))

        single_operands=split_expression[i].split(',')

        if split_expression[max(0,i-1)]=='=':
          if i-2<0:
            raise ValueError('Cannot evaluate subexpression '+' '.join(split_expression[max(0,i-1):min(i+1,len(split_expression))])+' starting with an \'=\'')
          main_operand=split_expression[i-2]
          temp_split_expression.pop()
          temp_split_expression.pop()
          skips=0

        else:
          if i+2>len(split_expression)-1:
            raise ValueError('Cannot evaluate subexpression '+' '.join(split_expression[max(0,i-1):min(i+1,len(split_expression))])+' ending with an \'=\'')
          main_operand=split_expression[i+2]
          skips=2

        temp_expression=list(['('])
        temp_expression+=' or '.join(['%s = %s'% (a,b) for (a,b) in zip(len(single_operands)*[main_operand],single_operands)]).split()
        temp_expression.append(')')
        temp_split_expression+=temp_expression
        continue

      temp_split_expression.append(split_expression[i])

    split_expression=temp_split_expression

    #handle ':' operator
    #replaces an expression like 'col_a=x:y' with '(col_a>=x and col_a<=y)'
        
    temp_split_expression=list()
    skips=0

    for i in range(len(split_expression)):
      if skips>0:
        skips-=1
        continue
      if ':' in split_expression[i]:
        if split_expression[max(0,i-1)] != '=' and split_expression[min(i+1,len(split_expression)-1)] != '=':
          raise ValueError('Can evaluate subexpression '+' '.join(split_expression[max(0,i-1):min(i+1,len(split_expression))])+' \':\' sign is only allowed in combination with \'=\'')
        if len(split_expression[i].split(':')) != 2:
          raise ValueError('Can operate \':\' operator only on 2 operands in subexpression '+' '.join(split_expression[max(0,i-1):min(i+1,len(split_expression))]))

        lhs=split_expression[i].split(':')[0]
        rhs=split_expression[i].split(':')[1]

        template_expression=['(','','<=','','and','','<=','',')']

        if split_expression[max(0,i-1)] == '=':
          if i-2<0:
            raise ValueError('Cannot evaluate subexpression '+' '.join(split_expression[max(0,i-1):min(i+1,len(split_expression))])+' starting with an \'=\'')
          temp_split_expression.pop()
          temp_split_expression.pop()
          template_expression[3]=split_expression[i-2]
          template_expression[5]=split_expression[i-2]
          skips=0

        else:
          if i+2>len(split_expression)-1:
            raise ValueError('Cannot evaluate subexpression '+' '.join(split_expression[max(0,i-1):min(i+1,len(split_expression))])+' ending with an \'=\'')
          template_expression[3]=split_expression[i+2]
          template_expression[5]=split_expression[i+2]
          skips=2 

        template_expression[1]=lhs
        template_expression[7]=rhs
        temp_split_expression+=template_expression
        continue

      temp_split_expression.append(split_expression[i])
  
    split_expression=temp_split_expression

    #the whole thing is now split to pieces, we need to cast the types of the operands
    #into the types of the columns, the operands are compared against.

    final_expression=list()
    subexpression=list()

    for item in split_expression:
      if item in self.valid_operators:
        if self.valid_operators[item]['type'] in ['boolean','left_bracket','right_bracket']:
          if len(subexpression)>0:
            #figure out the column type(s) and parse the operands
            final_expression+=self._ParseSubExpression(subexpression)
            subexpression=list()
          final_expression.append(item)
          continue
      subexpression.append(item)

    if len(subexpression)>0:
      final_expression+=self._ParseSubExpression(subexpression)

    return final_expression       
      

  def _ExpressionLexer(self, expression):

    #Reads token after token and searches for brackets and valid_operators
    #everything, that doesn't match the above is assumed to be an operand

    split_expression=list()

    actual_position=0
    eaten_stuff=''

    while True:

      if actual_position>=len(expression):
        if len(eaten_stuff)>0:
          split_expression.append(eaten_stuff)
        return split_expression

      token=expression[actual_position]

      if token.isspace():
        if len(eaten_stuff)>0:
          split_expression.append(eaten_stuff)
          eaten_stuff=''
        actual_position+=1
        continue

      #note, that there is no check for boolean operators. They need to be clearly separated by spaces
      #or brackets anyway, so they get added with the eaten stuff
      if token in self.valid_operators:
        if self.valid_operators[token]['type']=='left_bracket' or self.valid_operators[token]['type']=='right_bracket':
          if len(eaten_stuff)>0:
            split_expression.append(eaten_stuff)
            eaten_stuff=''
          split_expression.append(token)
          actual_position+=1
          continue

        if self.valid_operators[token]['type']=='arithmetic':
          if len(eaten_stuff)>0:
            split_expression.append(eaten_stuff)
            eaten_stuff=''
          split_expression.append(token)
          actual_position+=1
          continue

        if self.valid_operators[token]['type']=='comparison':
          if len(eaten_stuff)>0:
            split_expression.append(eaten_stuff)
            eaten_stuff=''
          if self.valid_operators[token]['extension']!=None:
            if actual_position+len(self.valid_operators[token]['extension'])<len(expression):
              if expression[actual_position:actual_position+len(self.valid_operators[token]['extension'])]==self.valid_operators[token]['extension']:
                split_expression.append(self.valid_operators[token]['extension'])
                actual_position+=len(self.valid_operators[token]['extension'])
                continue
          split_expression.append(token)
          actual_position+=1
          continue

      eaten_stuff+=token
      actual_position+=1
