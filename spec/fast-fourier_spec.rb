require_relative '../lib/FastFourier/fast-fourier.rb'

EPSILON = (2.0)**(-12.0)
module FFTHelper
  def FFTHelper.equal x, y
    if [x, y].any? {|a| a.is_a? Complex}
      (x - y).magnitude < EPSILON
    else
      (x - y) < EPSILON
    end
  end
end

include FastFourier

describe FastFourier do
  describe "root_of_unity" do
    2.times do
      den = 1 + rand(100)
      num = rand(den)

      ans = FastFourier.root_of_unity(den, num)
      specify "#{num}th root raised to #{num}th power should be idempotent" do
        #puts "#{ans**den}"
        expect((ans**den - 1).magnitude).to be < EPSILON
      end
    end
  end

  inputs = (0..6).map {|x| (0...(2**x)).map{rand(1000)} } + (5..23).step(8).map {|x| (0...x).map{rand(1000)}}

  describe "dft" do
    (0...inputs.length).each do |i|

      angle = (2 * Math::PI) / inputs[i].length
      npr = Complex(Math.cos(angle), Math.sin(angle))
      result = FastFourier.dft(inputs[i])
      expected_result = (0...(inputs[i].length)).map do |j|
        (0...(inputs[i].length)).map {|k| inputs[i][k] * npr ** (-1*j*k)}.inject(:+)
      end

      describe "for input of length #{inputs[i].length}" do 

        specify "the result should also be of length #{inputs[i].length}" do
          expect(result.length).to be == inputs[i].length
        end

        (0...(result.length)).each do |j|
          specify "index #{j} should match the DFT definition" do
            expect(result[j]).to be_within(EPSILON).of(expected_result[j])
          end
        end

        describe "the result should be invertable." do
          inv = FastFourier.dft(result, true)
          #puts "#{inversion} <=> #{inputs[i]}"

          (0...(result.length)).each do |j|
            specify "result index #{j} did not match input" do
              expect(inv[j]).to be_within(EPSILON).of(inputs[i][j])
            end
          end

          inv_c = FastFourier.dft(result, :conjugate)
          #puts "#{inv_c} <=> #{inputs[i]}"

          (0...(result.length)).each do |j|
            specify "conjugate result index #{j} did not match input" do
              expect(inv_c[j]).to be_within(EPSILON).of(inputs[i][j])
            end
          end

          inv_s = FastFourier.dft(result, :swap)
          #puts "#{inv_s} <=> #{inputs[i]}"

          (0...(result.length)).each do |j|
            specify "swap result index #{j} did not match input" do
              expect(inv_s[j]).to be_within(EPSILON).of(inputs[i][j])
            end
          end
        end
      end
    end
  end

  describe "discrete_involutary" do

    (0...inputs.length).each do |i|
      result = FastFourier.discrete_involutary(inputs[i])
      describe "for input of length #{inputs[i].length}" do 
        specify "the result should also be of length #{inputs[i].length}" do
          expect(result.length).to be == inputs[i].length
        end
        describe "the result should be invertable." do
          inv = FastFourier.discrete_involutary(result)
          (0...(result.length)).each do |j|
            specify "conjugate result index #{j} did not match input" do
              expect(inv[j]).to be_within(EPSILON).of(inputs[i][j])
            end
          end
        end
      end
    end
  end

end
