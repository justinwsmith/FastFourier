module FastFourier
  VERSION = "0.0.1"

  @@fast_fourier_roots = { Rational(0,1) => Complex(1, 0),
              Rational(1,4) => Complex(0, -1),  
              Rational(1,2) => Complex(-1, 0),  
              Rational(3,4) => Complex(0, 1),
              Rational(1,1) => Complex(1, 0)
            }

  def root_of_unity(den, num = 1)
    key = Rational(num, den) % 1
    if ans = @@fast_fourier_roots[key]
      ans
    else
      angle = -2*Math::PI * key
      @@fast_fourier_roots[key] = Complex(Math.cos(angle), Math.sin(angle))
    end
  end
  alias_method(:ff_root, :root_of_unity)

  def fast_fourier_extend(vals, inverse)
  	n = vals.length
  	if(n & n-1 > 0) # if n is not a power of 2
  		pow = Math.log2(n).ceil
  		pow2 = 2**pow
  		v = vals.clone
	  	if inverse
	  		(n...pow2).each{|i| v[i] = vals[pow2-i].conj }
	  	else
	  		(n...pow2).each{|i| v[i] = vals[i-n] }
	  	end
  		v
  	else
  		vals
  	end
  end

  def discrete_fourier_slow(vals, inverse=false)
  	npr = root_of_unity(vals.length)
  	(0...(vals.length)).map do |i| 
        ary = (0...(vals.length)).map do |j| 
        	vals[j] * npr ** (-1*i*j)
        end
        ary.inject(:+)
    end
  end

  def discrete_fourier(vals, inverse = false)
  	n = vals.length
		vals = fast_fourier_extend(vals, inverse)

    if inverse
      case inverse
      when :conjugate
        cooley_tukey_dft( vals.map { |x| x.conj } ).map {|x| x.conj / vals.length}
      when :swap
        cooley_tukey_dft( vals.map { |x| Complex(x.imag, x.real) } ).map {|x| Complex(x.imag, x.real) / vals.length}
      else
        cooley_tukey_dft([vals[0]] + vals[1..-1].reverse).map{|x| x / vals.length}
      end
    else
      cooley_tukey_dft(vals)
    end
  end
  alias_method(:dft, :discrete_fourier)

  def discrete_involutary(vals)
  	n = vals.length
		vals = fast_fourier_extend(vals, false)

    result = cooley_tukey_dft( vals.map { |x| x.conj } ).map {|x| x / Math.sqrt(vals.length) }

    result[0...n]    
  end

  def cooley_tukey_dft(vals)
    if vals.size == 1
      vals
    else
      n = vals.length
      vals_even = cooley_tukey_dft((0...n).step(2).map{|i| vals[i]})
      vals_odd = cooley_tukey_dft((1...n).step(2).map{|i| vals[i]})
      ret_vals = Array.new(vals.size)
      #npr = NPRs[n] || Math::E**((-2 * Math::PI * Complex(0, 1)) / n)
      (0...(n/2)).each do |i|
        ei = vals_even[i]
        oi = vals_odd[i]
        npri = root_of_unity(n, i)
        oink = npri * oi
        ret_vals[i] = (ei + oink)
        ret_vals[i+n/2] = (ei - oink)
      end
      ret_vals
    end
  end

end
