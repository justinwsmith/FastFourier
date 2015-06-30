=begin
The MIT License (MIT)

Copyright (c) 2014 Justin W, Smith

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
=end

module FastFourier

  @@fast_fourier_roots = {
    Rational(0,1) => Complex(1, 0),
    Rational(1,4) => Complex(0, 1),
    Rational(1,2) => Complex(-1, 0),
    Rational(3,4) => Complex(0, -1),
    Rational(1,1) => Complex(1, 0)
  }

  def FastFourier.discrete_fourier(vals, inverse = false)
    n = vals.length
    if(n & (n-1) == 0)
      dft =  :cooley_tukey_dft
    else
      dft = :discrete_fourier_slow
    end

    if inverse
      case inverse
        when :conjugate
          send( dft, vals.map { |x| x.conj } ).map {|x| x.conj / vals.length}
        when :swap
          send( dft,  vals.map { |x| Complex(x.imag, x.real) } ).map {|x| Complex(x.imag, x.real) / vals.length}
        else
          send( dft, [vals[0]] + vals[1..-1].reverse).map{|x| x / vals.length}
      end
    else
      send( dft, vals)
    end
  end

  class <<self
    alias_method(:dft, :discrete_fourier)
  end

  # This is a self-inverse variation of the discrete fourier transform
  def FastFourier.discrete_involutary(vals)
  	n = vals.length
    if(n % 2 == 0)
      dft =  :cooley_tukey_dft
    else
      dft = :discrete_fourier_slow
    end

    result = send( dft, vals.map { |x| x.conj } ).map {|x| x / Math.sqrt(vals.length) }

    result[0...n]    
  end

private

  def FastFourier.root_of_unity(den, num = 1)
    key = Rational(num, den) % 1
    if ans = @@fast_fourier_roots[key]
      ans
    else
      angle = 2*Math::PI * key
      @@fast_fourier_roots[key] = Complex(Math.cos(angle), Math.sin(angle))
    end
  end

  def FastFourier.cooley_tukey_dft(vals)
    if vals.size == 1
      vals
    else
      n = vals.length
      if(n % 2 == 0) # if n is not a power of 2
        dft =  :cooley_tukey_dft
      else
        dft = :discrete_fourier_slow
      end
      vals_even = send( dft, (0...n).step(2).map{|i| vals[i]})
      vals_odd = send( dft, (1...n).step(2).map{|i| vals[i]})

      ret_vals = Array.new(n)

      (0...(n/2)).each do |i|
        ei = vals_even[i]
        oi = vals_odd[i]
        npri = root_of_unity(n, i).conj
        oink = npri * oi
        ret_vals[i] = (ei + oink)
        ret_vals[i+n/2] = (ei - oink)
      end
      ret_vals
    end
  end

  def FastFourier.discrete_fourier_slow(vals, inverse=false)
    npr = root_of_unity(vals.length)
    (0...(vals.length)).map do |i|
      ary = (0...(vals.length)).map do |j|
        vals[j] * npr ** (-1*i*j)
      end
      ary.inject(:+)
    end
  end

end
