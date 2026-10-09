# Numerical Lab

เว็บไซต์เรียนรู้ Numerical Methods ภาษาไทย เน้น visualization ใช้ HTML, CSS และ JavaScript แบบ OOP โดยไม่ต้องติดตั้ง dependencies

เปิด `dist/index.html` ใน browser ได้ทันที หรือให้ HTTP server serve โฟลเดอร์ `dist`

## วิธีคำนวณ

- Bisection: ตรวจว่าปลายช่วงเปลี่ยนเครื่องหมาย หยุดเมื่อ residual หรือครึ่งความกว้างช่วงถึง tolerance
- Newton–Raphson: คำนวณอนุพันธ์ด้วย central difference หยุดตาม residual
- Secant: ใช้สองจุดประมาณราก หยุดตาม residual
- Composite Trapezoidal และ Simpson's 1/3: แสดงพื้นที่สะสม Simpson ต้องใช้จำนวนช่วงเป็นเลขคู่
- Euler: ประมาณคำตอบ dy/dx = f(x,y) จากค่าเริ่มต้น
- Conjugate Gradient: แก้ระบบ Ax = b ขนาด 2 × 2 โดยตรวจว่า A สมมาตรและ positive definite แสดงเส้นชั้นความสูงของ quadratic objective เส้นทางค้นหา และ residual แบบ log scale พร้อม α และ β ในตาราง

สำหรับปริพันธ์และ Euler ตัวเลขผลต่าง n กับ 2n เป็นการเปรียบเทียบความละเอียด ไม่ใช่ขอบเขต error ที่รับประกัน สำหรับ Bisection ต้องใช้ฟังก์ชันต่อเนื่องในช่วงเพื่อให้การรับประกันการหารากมีผล

## โครงสร้าง OOP

`Expression` แปลงสมการเป็น syntax tree โดยไม่ใช้ eval หรือ Function

`NumericalMethod` เป็นคลาสแม่ของ `Bisection`, `NewtonRaphson`, `Secant`, `Quadrature`, `EulerMethod` แยกการคำนวณออกจาก DOM

`PlotRenderer` วาดกราฟบน Canvas และปรับตามขนาดหน้าจอ

`NumericalLabApp` จัดการฟอร์ม ตาราง animation และ CSV export

`ConjugateGradient` ใน `dist/conjugate-gradient.js` สืบทอดจาก `NumericalMethod` และแยกการคำนวณออกจาก `ConjugateGradientPlot` กับ `ExtendedNumericalLabApp` ใช้ residual จริง ‖b − Ax‖₂ เป็นเกณฑ์หยุด ค่าเริ่มต้นที่ผ่าน tolerance แล้วจะแสดง 0 รอบ

รองรับ + - * / ^, วงเล็บ, x, y (Euler), pi, e, sin, cos, tan, exp, log, sqrt, abs โดยมุมเป็นเรเดียน ฟอนต์จาก Google Fonts มี system fallback เมื่อ offline

## ตรวจสอบ

ใช้ `node --test tests/numerical.test.cjs tests/conjugate-gradient.test.cjs` ตรวจความถูกต้องของตัวอย่างและเงื่อนไขข้อผิดพลาด รวมทั้งการเป็น A-conjugate และกรณีระบบไม่เป็น SPD
